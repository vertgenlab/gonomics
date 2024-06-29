package fileio

import (
	"bufio"
	"errors"
	"fmt"
	"io"
	"unicode/utf8"
)

// Read implements io.Reader, reading data into b from the file or gzip stream.
func (br *ByteReader) Read(b []byte) (n int, err error) {
    // GZIP reading with pooling
    if br.internalGzip != nil {
        return br.internalGzip.Read(b) // Directly read into the output buffer
    } else {
		// Non-GZIP reading with pooling
		return br.Reader.Read(b)
	}
    
}


// Close closes the internal File and gzip Reader if present.
func (br *ByteReader) Close() error {
	if br.File == nil {
		return errors.New("no file found")
	}
	var err error
	if br.internalGzip != nil {
		err = br.internalGzip.Close()
	}
	if fileErr := br.File.Close(); fileErr != nil {
		if err != nil {
			return fmt.Errorf("multiple errors: gzip: %w, file: %w", err, fileErr)
		}
		err = fmt.Errorf("closing file: %w", fileErr) // Wrap the file error for clarity
	}
	return err
}

// getBuffer retrieves a xbuffer from the pool.
func (b *ByteReader) getBuffer() []byte {
	buf := b.bufPool.Get().([]byte)
	buf = buf[:0]
	return buf
}

// getBuffer retrieves a xbuffer from the pool.
func (b *ByteWriter) getBuffer() []byte {
	return b.bufPool.Get().([]byte)
}


// putBuffer returns a buffer to the pool.
func (b *ByteReader) putBuffer(buf []byte) {
	b.bufPool.Put(buf)
}

func NewWriterAutoFlush(w io.Writer, size int, flushAt float32) *ByteWriter {
	if flushAt < 0 && flushAt != -1 || flushAt > 1 {
		panic(fmt.Sprintf("flushAt should be -1.0 or a fraction 0.0 <= flushAt <= 1.0, got %f", flushAt))
	}
	// Is it already a ByteWriter?
	_, ok := w.(*ByteWriter)
	if !ok {
		_, ok = w.(*ByteWriter)
	}
	if ok {
		panic("Will not auto flush on top of a buffered writer")
	}

	b := NewByteWriterSize(w, size)
	if flushAt != -1 {
		b.flushAt = int(flushAt * float32(size))
	}
	return b
}

// Reset discards any unflushed buffered data, clears any error, and
// resets b to write its output to w. If w.err is nil, in order to ensure
// no partial writes end up in w, it waits until any chunked write and/or
// flush complete before the output is redirected.
func (b *ByteWriter) Reset(w *bufio.Writer) {
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

// Flush writes any buffered data to the underlying io.ByteWriter. As long as the
// buffer has enough available space, writes can proceed concurrently.
func (b *ByteWriter) Flush() error {
	b.mtx.Lock()
	b.waitUntilNoChunkedWrite()

	err := b.flush(0)

	b.mtx.Unlock()
	return err
}

// Flushing implementation: need is the minimum availabla buffer space expected
// by the caller.
func (b *ByteWriter) flush(need int) error {
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
func (b *ByteWriter) maybeAutoFlush() {
	if b.n >= b.flushAt && b.nFlush == 0 {
		go b.Flush()
	}
}

// Resets the inChunkedWriteMode flag and wakes a goroutine waiting to write.
func (b *ByteWriter) endChunkedWrite() {
	b.inChunkedWriteMode = false
	// Wake all goroutines waiting to write
	b.noChunkedWrite.Signal()
}

// Waits for any in-progress chunked write to complete. This is safe to do even
// in the face of errors, as a chunked write will reset the chunked write mode
// flag and signal moChunkedWrite regardless of error.
func (b *ByteWriter) waitUntilNoChunkedWrite() {
	for b.inChunkedWriteMode {
		b.noChunkedWrite.Wait()
	}
	b.noChunkedWrite.Signal()
}

// Available returns how many bytes are unused in the buffer.
func (b *ByteWriter) Available() int {
	b.mtx.Lock()
	res := b.available()
	b.mtx.Unlock()
	return res
}
func (b *ByteWriter) available() int { return len(b.buf) - b.n }

// Buffered returns the number of bytes that have been written into the current buffer.
func (b *ByteWriter) Buffered() int {
	b.mtx.Lock()
	res := b.buffered()
	b.mtx.Unlock()
	return res
}
func (b *ByteWriter) buffered() int { return b.n }

// Write writes the contents of p into the buffer.
// It returns the number of bytes written.
// If nn < len(p), it also returns an error explaining
// why the write is short.
func (b *ByteWriter) Write(p []byte) (nn int, err error) {
	b.mtx.Lock()
	defer b.mtx.Unlock()
	b.waitUntilNoChunkedWrite()
	
    if b.internalGzip != nil {
        return b.internalGzip.Write(p)
    }

    return b.wr.Write(p)
}


// WriteByte writes a single byte.
func (b *ByteWriter) WriteByte(c byte) error {
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
func (b *ByteWriter) WriteRune(r rune) (size int, err error) {
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
func (b *ByteWriter) WriteString(s string) (int, error) {
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