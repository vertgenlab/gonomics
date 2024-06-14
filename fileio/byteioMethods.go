package fileio

import (
	"errors"
	"io"
	"log"
)

// Read implements io.Reader, reading data into b from the file or gzip stream.
func (reader *ByteReader) Read(b []byte) (n int, err error) {
	if reader.internalGzip != nil {
		return reader.internalGzip.Read(b)
	}
	return reader.Reader.Read(b)
}

// Close closes the internal File and gzip Reader if present.
func (br *ByteReader) Close() error {
	var gzErr, fileErr error
	if br.internalGzip != nil {
		gzErr = br.internalGzip.Close()
	}
	if br.File != nil {
		fileErr = br.File.Close()
	} else {
		return errors.New("no file found")
	}

	switch { // Handle error returns. Priority is gzErr > fileErr
	case gzErr != nil:
		return gzErr

	case fileErr != nil:
		log.Println("WARNING: attempted to close file, but file already closed")
		return nil

	default:
		return nil
	}
}

// Write writes the contents of p into the buffer.
func (b *ByteWriter) Write(p []byte) (nn int, err error) {
	if b.internalGzip != nil {
		return b.internalGzip.Write(p)
	}
	b.mtx.Lock()
	defer b.mtx.Unlock()
	if b.closed {
		return 0, io.ErrClosedPipe
	}
	b.waitUntilNoChunkedWrite()

	for len(p) > len(b.buf)-b.n && b.err == nil {
		if !b.inChunkedWriteMode {
			b.inChunkedWriteMode = true
			defer b.endChunkedWrite()
		}

		var n int
		if b.n == 0 {
			n, b.err = b.Writer.Write(p)
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
func (b *ByteWriter) WriteByte(c byte) error {
	b.mtx.Lock()
	defer b.mtx.Unlock()
	if b.closed {
		return io.ErrClosedPipe
	}
	b.waitUntilNoChunkedWrite()

	if b.err != nil {
		return b.err
	}

	if len(b.buf)-b.n <= 0 && b.flush(1) != nil {
		return b.err
	}
	b.buf[b.n] = c
	b.n++
	b.maybeAutoFlush()
	return nil
}

// WriteString writes a string.
func (b *ByteWriter) WriteString(s string) (int, error) {
	if b.internalGzip != nil {
		return b.internalGzip.Write([]byte(s))
	}
	nn := 0
	b.mtx.Lock()
	defer b.mtx.Unlock()
	b.waitUntilNoChunkedWrite()

	if b.err != nil {
		return 0, b.err
	}

	for len(s) > len(b.buf)-b.n {
		n := copy(b.buf[b.n:], s)
		b.n += n
		nn += n
		s = s[n:]
		if !b.inChunkedWriteMode {
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

// Close flushes any buffered data and releases resources.
func (b *ByteWriter) Close() error {
	if b.internalGzip != nil {
		return b.internalGzip.Close()
	}
	b.mtx.Lock()
	defer b.mtx.Unlock()

	if b.closed {
		return io.ErrClosedPipe
	}

	// Wait for any in-progress chunked writes or flushes to complete
	for b.err == nil && (b.nFlush != 0 || b.inChunkedWriteMode) {
		b.notFlushing.Wait()
	}

	// Perform final flush
	err := b.flush(0)

	// Return buffer to pool
	b.putBuffer(b.buf)
	b.closed = true

	return err
}

// Flush writes any buffered data to the underlying io.Writer.
func (b *ByteWriter) Flush() error {
	b.mtx.Lock()
	defer b.mtx.Unlock()
	b.waitUntilNoChunkedWrite()

	return b.flush(0)
}

// flush writes the buffered data to the underlying io.Writer.
func (b *ByteWriter) flush(need int) error {
	iAmChunked := b.inChunkedWriteMode

	callerNeedMet := func() bool {
		return b.n == 0 || (need > 0 && len(b.buf)-b.n >= need)
	}

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

	for b.err == nil && !callerNeedMet() {
		if waitToFlush() {
			continue
		}

		b.nFlush = b.n

		mtxReleased := false
		if !iAmChunked && b.nFlush != len(b.buf) {
			b.mtx.Unlock()
			mtxReleased = true
		}

		n, err := b.Writer.Write(b.buf[:b.nFlush])
		if n < b.nFlush && err == nil {
			err = io.ErrShortWrite
		}

		if mtxReleased {
			b.mtx.Lock()
		}

		if n > 0 && n < b.n {
			copy(b.buf[:b.n-n], b.buf[n:b.n])
		}
		b.n -= n
		b.err = err
		b.nFlush = 0

		if mtxReleased && b.inChunkedWriteMode {
			b.chunkedWriter.Signal()
			b.waitUntilNoChunkedWrite()
		}

		if need == 0 {
			break
		}
	}

	b.notFlushing.Signal()
	return b.err
}

// maybeAutoFlush triggers an async Flush() if more than flushAt bytes are used and no Flush() call is in progress.
func (b *ByteWriter) maybeAutoFlush() {
	if b.n >= b.flushAt && b.nFlush == 0 {
		go b.Flush()
	}
}

// Reset discards any unflushed buffered data, clears any error, and resets b to write its output to w.
func (b *ByteWriter) Reset(w io.Writer) {
	b.mtx.Lock()
	defer b.mtx.Unlock()

	for b.err == nil && (b.nFlush != 0 || b.inChunkedWriteMode) {
		b.notFlushing.Wait()
	}

	b.err = nil
	b.n = 0
	b.Writer = w
	b.nFlush = 0

	b.notFlushing.Signal()
}

// Buffered returns the number of bytes that have been written into the current buffer.
func (b *ByteWriter) Buffered() int {
	b.mtx.Lock()
	defer b.mtx.Unlock()
	return b.n
}

// getBuffer retrieves a buffer from the pool.
func (b *ByteWriter) getBuffer() []byte {
	return b.bufPool.Get().([]byte)
}

// putBuffer returns a buffer to the pool.
func (b *ByteWriter) putBuffer(buf []byte) {
	b.bufPool.Put(buf)
}

// endChunkedWrite resets the inChunkedWriteMode flag and wakes a goroutine waiting to write.
func (b *ByteWriter) endChunkedWrite() {
	b.inChunkedWriteMode = false
	b.noChunkedWrite.Signal()
}

// waitUntilNoChunkedWrite waits for any in-progress chunked write to complete.
func (b *ByteWriter) waitUntilNoChunkedWrite() {
	for b.inChunkedWriteMode {
		b.noChunkedWrite.Wait()
	}
	b.noChunkedWrite.Signal()
}
