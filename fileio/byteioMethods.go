package fileio

import (
	"errors"
	"fmt"
	"io"
)

// Read implements io.Reader, reading data into b from the file or gzip stream.
func (br *ByteReader) Read(b []byte) (n int, err error) {
	// Buffer empty, get a new one from the pool
	if br.buf != nil {
		br.putBuffer(br.buf) // Return the old buffer to the pool
	}
	br.buf = br.getBuffer()
	defer br.putBuffer(br.buf)
	// GZIP reading with pooling
	if br.internalGzip != nil {
		n, err = br.internalGzip.Read(br.buf)
		if err != nil {
			return n, err
		}

		// Copy decompressed data from tempBuf into the provided buffer 'b'
		n = copy(b, br.buf[:n])
		return n, nil
	}
	// Non-GZIP reading with pooling
	for len(b) > 0 {
		if len(br.buf) > 0 {
			// Data available in existing buffer
			copied := copy(b, br.buf)
			br.buf = br.buf[copied:]
			b = b[copied:]
			n += copied
		} else {

			nr, err := br.Reader.Read(br.buf)
			if nr > 0 {
				copied := copy(b, br.buf)
				br.buf = br.buf[copied:]
				b = b[copied:]
				n += copied
			}
			if err != nil {
				return n, err
			}
		}
	}
	return n, nil
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

// getBuffer retrieves a xbuffer from the pool.
func (b *ByteReader) getBuffer() []byte {
	return b.bufPool.Get().([]byte)
}

// getBuffer retrieves a buffer from the pool.
func (b *ByteWriter) getBuffer() []byte {
	return b.bufPool.Get().([]byte)
}

// putBuffer returns a buffer to the pool.
func (b *ByteReader) putBuffer(buf []byte) {
	b.bufPool.Put(buf)
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
