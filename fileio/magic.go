package fileio

import (
	"bytes"
	"io"
	"os"

	"github.com/vertgenlab/gonomics/exception"
)

// magic number that identifies a gzip file
var magicGzip = []byte{0x1f, 0x8b}

// IsGzip checks if the magic gzip number is
// present at the start of the input io.ReadSeeker.
// After check for the magic number, IsGzip
// seeks back to the initial position.
func IsGzip(r io.ReadSeeker) bool {
	return checkMagic(r, magicGzip)
}

// checkMagic returns true if the input io.ReadSeeker
// begins with the bytes in b. Seeks back to initial
// position after checking for b.
func checkMagic(r io.ReadSeeker, b []byte) bool {
	// get current pos
	initialPos, err := r.Seek(0, io.SeekCurrent)
	exception.PanicOnErr(err)

	// seek to start
	_, err = r.Seek(0, io.SeekStart)
	exception.PanicOnErr(err)

	readBytes := make([]byte, len(b))
	r.Read(readBytes)

	// seek back to initialPos
	_, err = r.Seek(initialPos, io.SeekStart)
	exception.PanicOnErr(err)

	// check readBytes against b
	return bytes.Equal(readBytes, b)
}

// newStdinMagicReader creates a reader with an internal buffer of len(magic) that
// serves the read magic bytes to a function calling Read, then serves the rest of
// os.Stdin without passing the read bytes through the internal buffer.
//
// This basically allows the function to peek at the first bytes of the stdin stream
// to see if they match the input magic bytes.
func newStdinMagicReader(magic []byte) (reader *stdinMagicReader, hasMagic bool) {
	// create reader
	reader = &stdinMagicReader{
		in:             os.Stdin,
		readMagicBytes: make([]byte, len(magic)),
	}

	// read len(magic) bytes into internal buffer
	// and truncate internal buffer if necessary
	n, _ := reader.in.Read(reader.readMagicBytes)
	reader.readMagicBytes = reader.readMagicBytes[:n]

	// check if internal buffer matches input magic bytes
	hasMagic = true
	if n != len(magic) {
		hasMagic = false
	} else {
		hasMagic = bytes.Equal(magic, reader.readMagicBytes)
	}

	return
}

// stdinMagicReader is designed to have a small internal buffer which
// is filled at creation and served on the first call to Read with
// further calls to Read reading from os.Stdin directly.
type stdinMagicReader struct {
	in             *os.File
	readMagicBytes []byte
}

// Read will serve the bytes in readMagicBytes before serving the rest
// of os.Stdin directly.
func (r *stdinMagicReader) Read(b []byte) (n int, err error) {
	if r.readMagicBytes == nil {
		return r.in.Read(b)
	}

	mLen := len(r.readMagicBytes)
	switch {
	case len(b) > mLen:
		copy(b[0:mLen], r.readMagicBytes)
		n, err = r.in.Read(b[mLen:])
		n += mLen
		r.readMagicBytes = nil
		return

	case len(b) < mLen:
		copy(b, r.readMagicBytes)
		n = len(b)
		r.readMagicBytes = r.readMagicBytes[len(b):]
		return

	default: // len(b) == mLen
		copy(b, r.readMagicBytes)
		r.readMagicBytes = nil
		return
	}
}
