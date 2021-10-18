package fileio

import (
	"errors"
	"log"
)

// Close the receiving EasyReader.
func (er *EasyReader) Close() error {
	var gzErr, fileErr error
	if er.internalGzip != nil {
		gzErr = er.internalGzip.Close()
	}
	if er.File != nil {
		fileErr = er.File.Close()
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

// Read retrieves n bytes from the receiving EasyReader.
func (er *EasyReader) Read(p []byte) (n int, err error) {
	if er.internalGzip != nil {
		return er.internalGzip.Read(p)
	} else {
		return er.BuffReader.Read(p)
	}
}

// Close the receiving EasyWriter.
func (ew *EasyWriter) Close() error {
	var gzErr, bufErr, fileErr error
	if ew.internalGzip != nil {
		gzErr = ew.internalGzip.Close() // Serious write errors possible.
	}
	if ew.internalBuff != nil {
		bufErr = ew.internalBuff.Flush() // Serious write errors possible.
	}
	if ew.File != nil {
		fileErr = ew.File.Close() // The only possible err is that the file has already been closed.
	} else {
		return errors.New("no open file")
	}

	switch { // Handle error returns. Priority is gzErr > bufErr > fileErr
	case gzErr != nil:
		return gzErr

	case bufErr != nil:
		return bufErr

	case fileErr != nil:
		log.Println("WARNING: attempted to close file, but file already closed")
		return nil

	default:
		return nil
	}
}

// Write bytes to the receiving EasyWriter.
func (ew *EasyWriter) Write(p []byte) (n int, err error) {
	if ew.internalGzip != nil {
		return ew.internalGzip.Write(p)
	} else {
		return ew.internalBuff.Write(p)
	}
}

// Peek retrieves the next n bytes of the receiving EasyReader without advancing the reader.
func (er *EasyReader) Peek(n int) ([]byte, error) {
	return er.BuffReader.Peek(n)
}
