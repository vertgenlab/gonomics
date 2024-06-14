package fileio

import (
	"bufio"
	"fmt"
	"io"
	"sync"
	"unicode/utf8"
)

const maxConsecutiveEmptyReads = 100

// Writer implements highly concurrent buffering for an io.Writer object.
// In particular, writes will not block while a Flush() call is in progress as

type Writer struct {
	// Protects all internal state
	mtx *sync.Mutex

	err error
	buf []byte
	n   int
	wr  io.Writer

	// Whether a chunked write is in flight and Writer is in "chunked write mode"
	inChunkedWriteMode bool
	// Condition variable for all goroutines waiting on a chunked write
	noChunkedWrite *sync.Cond

	// Fields below are only used by flush() (and marginally by Reset()) to
	// provide correct serialization of writes.

	// Number of bytes currently being flushed from the start of buf. Only
	// non-zero while a flush is in progress. Always <= n.
	nFlush int
	// Condition variable for all goroutines waiting to flush()
	notFlushing *sync.Cond
	// Condition variable that the only chunked writer is waiting on to flush()
	chunkedWriter *sync.Cond

	// Automatically flush when n > flushAt
	flushAt int
}

// NewWriterSize returns a new Writer whose buffer has at least the specified
// size. If the argument io.Writer is already a Writer with large enough size,
// it returns the underlying Writer.
func NewWriterSize(w io.Writer, size int) *Writer {
	// Is it already a Writer?
	b, ok := w.(*Writer)
	if ok && len(b.buf) >= size {
		return b
	}
	if size <= 0 {
		size = defaultBufSize
	}
	if size < utf8.UTFMax {
		size = utf8.UTFMax
	}
	m := new(sync.Mutex)
	return &Writer{
		mtx:            m,
		buf:            make([]byte, size),
		wr:             w,
		noChunkedWrite: sync.NewCond(m),
		notFlushing:    sync.NewCond(m),
		chunkedWriter:  sync.NewCond(m),
		flushAt:        2 * size,
	}
}

// NewWriter returns a new Writer whose buffer has the default size.
func NewWriter(w io.Writer) *Writer {
	return NewWriterSize(w, defaultBufSize)
}

// NewWriterAutoFlush returns a new Writer whose buffer has at least the
// specified size and that will automatically trigger an asynchronous flush when
// the given buffer fraction is filled (e.g. 0.75 will flush when the buffer is
// 75% full). Panics if the argument io.Writer is already a Writer or
// bufio.Writer.
func NewWriterAutoFlush(w io.Writer, size int, flushAt float32) *Writer {
	if flushAt < 0 && flushAt != -1 || flushAt > 1 {
		panic(fmt.Sprintf("flushAt should be -1.0 or a fraction 0.0 <= flushAt <= 1.0, got %f", flushAt))
	}
	// Is it already a Writer?
	_, ok := w.(*Writer)
	if !ok {
		_, ok = w.(*bufio.Writer)
	}
	if ok {
		panic("Will not auto flush on top of a buffered writer")
	}

	b := NewWriterSize(w, size)
	if flushAt != -1 {
		b.flushAt = int(flushAt * float32(size))
	}
	return b
}

// Reset discards any unflushed buffered data, clears any error, and
// resets b to write its output to w. If w.err is nil, in order to ensure
// no partial writes end up in w, it waits until any chunked write and/or
// flush complete before the output is redirected.
func (b *Writer) Reset(w io.Writer) {
	b.mtx.Lock()

	// Wait until there are no chunked writes or flushes in progress
	for b.err == nil && (b.nFlush != 0 || b.inChunkedWriteMode) {
		b.notFlushing.Wait()
	}

	b.err = nil
	b.n = 0
	b.wr = w
	b.nFlush = 0
	b.mtx.Unlock()

	// (Potentially) wake one other goroutine waiting on flush()
	b.notFlushing.Signal()
}

// Flush writes any buffered data to the underlying io.Writer. As long as the
// buffer has enough available space, writes can proceed concurrently.
func (b *Writer) Flush() error {
	b.mtx.Lock()
	b.waitUntilNoChunkedWrite()

	err := b.flush(0)

	b.mtx.Unlock()
	return err
}

// Flushing implementation: need is the minimum availabla buffer space expected
// by the caller.
func (b *Writer) flush(need int) error {
	iAmChunked := b.inChunkedWriteMode

	callerNeedMet := func() bool {
		return b.n == 0 || (need > 0 && b.available() >= need)
	}

	// Waits for any in-progress flush() to complete, returns true iff mtx was
	// released and reaxquired in the process (i.e. if any Wait() call was made).
	waitToFlush := func() bool {
		waited := false
		if iAmChunked {
			for b.nFlush != 0 {
				b.chunkedWriter.Wait()
				waited = true
			}
		} else {
			for b.nFlush != 0 || b.inChunkedWriteMode {
				b.notFlushing.Wait()
				waited = true
			}
		}
		return waited
	}

	// Loop as long as no error and caller need is not yet met
	for b.err == nil && !callerNeedMet() {
		if waitToFlush() {
			continue
		}

		b.nFlush = b.n

		mtxReleased := false
		if !iAmChunked && b.nFlush != len(b.buf) {
			// Release the mutex to allow concurrent writes
			b.mtx.Unlock()
			mtxReleased = true
		}

		// Actually flush the first nFlush bytes
		n, err := b.wr.Write(b.buf[0:b.nFlush])
		if n < b.nFlush && err == nil {
			err = io.ErrShortWrite
		}

		if mtxReleased {
			// Grab back the mutex once the potentially blocking I/O call is done
			b.mtx.Lock()
		}

		if n > 0 && n < b.n {
			copy(b.buf[0:b.n-n], b.buf[n:b.n])
		}
		b.n -= n
		b.err = err
		b.nFlush = 0

		if mtxReleased && b.inChunkedWriteMode {
			// We are not the chunked writer but one exists: wake it and wait until
			// it is done (else our caller will enter into a race condition with it).
			b.chunkedWriter.Signal()
			b.waitUntilNoChunkedWrite()
		}

		if need == 0 {
			// Flush() call, no specific buffer space requirement, all done
			break
		}
	}

	// Always (potentially) wake a goroutine waiting to flush()
	b.notFlushing.Signal()
	return b.err
}

// Triggers an async Flush() if more than flushAt bytes are used and no Flush()
// call is already in progress.
func (b *Writer) maybeAutoFlush() {
	if b.n >= b.flushAt && b.nFlush == 0 {
		go b.Flush()
	}
}

// Resets the inChunkedWriteMode flag and wakes a goroutine waiting to write.
func (b *Writer) endChunkedWrite() {
	b.inChunkedWriteMode = false
	// Wake all goroutines waiting to write
	b.noChunkedWrite.Signal()
}

// Waits for any in-progress chunked write to complete. This is safe to do even
// in the face of errors, as a chunked write will reset the chunked write mode
// flag and signal moChunkedWrite regardless of error.
func (b *Writer) waitUntilNoChunkedWrite() {
	for b.inChunkedWriteMode {
		b.noChunkedWrite.Wait()
	}
	b.noChunkedWrite.Signal()
}

// Available returns how many bytes are unused in the buffer.
func (b *Writer) Available() int {
	b.mtx.Lock()
	res := b.available()
	b.mtx.Unlock()
	return res
}
func (b *Writer) available() int { return len(b.buf) - b.n }

// Buffered returns the number of bytes that have been written into the current buffer.
func (b *Writer) Buffered() int {
	b.mtx.Lock()
	res := b.buffered()
	b.mtx.Unlock()
	return res
}
func (b *Writer) buffered() int { return b.n }

// Write writes the contents of p into the buffer.
// It returns the number of bytes written.
// If nn < len(p), it also returns an error explaining
// why the write is short.
func (b *Writer) Write(p []byte) (nn int, err error) {
	b.mtx.Lock()
	defer b.mtx.Unlock()
	b.waitUntilNoChunkedWrite()

	for len(p) > b.available() && b.err == nil {
		if !b.inChunkedWriteMode {
			// Enter chunked write mode
			b.inChunkedWriteMode = true
			defer b.endChunkedWrite()
		}

		var n int
		if b.buffered() == 0 {
			// Large write, empty buffer.
			// Write directly from p to avoid copy.
			n, b.err = b.wr.Write(p)
		} else {
			n = copy(b.buf[b.n:], p)
			b.n += n
			b.flush(1)
		}
		nn += n
		p = p[n:]
	}
	if b.err != nil {
		return nn, b.err
	}
	n := copy(b.buf[b.n:], p)
	b.n += n
	nn += n
	b.maybeAutoFlush()
	return nn, nil
}

// WriteByte writes a single byte.
func (b *Writer) WriteByte(c byte) error {
	b.mtx.Lock()
	defer b.mtx.Unlock()
	b.waitUntilNoChunkedWrite()

	if b.err != nil {
		return b.err
	}

	if b.available() <= 0 && b.flush(1) != nil {
		return b.err
	}
	b.buf[b.n] = c
	b.n++
	b.maybeAutoFlush()
	return nil
}

// WriteRune writes a single Unicode code point, returning
// the number of bytes written and any error.
func (b *Writer) WriteRune(r rune) (size int, err error) {
	if r < utf8.RuneSelf {
		err = b.WriteByte(byte(r))
		if err != nil {
			return 0, err
		}
		return 1, nil
	}

	var encoded [4]byte
	size = utf8.EncodeRune(encoded[:], r)
	return b.Write(encoded[:size])
}

// WriteString writes a string.
// It returns the number of bytes written.
// If the count is less than len(s), it also returns an error explaining
// why the write is short.
func (b *Writer) WriteString(s string) (int, error) {
	nn := 0
	b.mtx.Lock()
	defer b.mtx.Unlock()
	b.waitUntilNoChunkedWrite()

	if b.err != nil {
		return 0, b.err
	}

	for len(s) > b.available() {
		n := copy(b.buf[b.n:], s)
		b.n += n
		nn += n
		s = s[n:]
		if !b.inChunkedWriteMode {
			// Enter chunked write mode
			b.inChunkedWriteMode = true
			defer b.endChunkedWrite()
		}
		b.flush(1)
		if b.err != nil {
			return nn, b.err
		}
	}
	n := copy(b.buf[b.n:], s)
	b.n += n
	nn += n
	b.maybeAutoFlush()
	return nn, nil
}

// ReadFrom implements io.ReaderFrom.
func (b *Writer) ReadFrom(r io.Reader) (n int64, err error) {
	b.mtx.Lock()
	defer b.mtx.Unlock()
	b.waitUntilNoChunkedWrite()

	if b.buffered() == 0 {
		if w, ok := b.wr.(io.ReaderFrom); ok {
			return w.ReadFrom(r)
		}
	}

	var m int
	for {
		if b.available() == 0 {
			if err1 := b.flush(1); err1 != nil {
				return n, err1
			}
		}
		nr := 0
		for nr < maxConsecutiveEmptyReads {
			m, err = r.Read(b.buf[b.n:])
			if m != 0 || err != nil {
				break
			}
			nr++
		}
		if nr == maxConsecutiveEmptyReads {
			b.maybeAutoFlush()
			return n, io.ErrNoProgress
		}
		b.n += m
		n += int64(m)
		if err != nil {
			break
		}
		if !b.inChunkedWriteMode {
			// Enter chunked write mode
			b.inChunkedWriteMode = true
			defer b.endChunkedWrite()
		}
	}
	if err == io.EOF {
		err = nil
	}
	if err != nil {
		b.maybeAutoFlush()
	}
	return n, err
}
