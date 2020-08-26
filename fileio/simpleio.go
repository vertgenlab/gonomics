package fileio

import(
	"strings"
	"os"
	"bufio"
	"io"
	"log"
	"compress/gzip"
	"github.com/vertgenlab/gonomics/common"
)

type SimpleReader struct {
	*bufio.Reader
	File *os.File
}

func (reader *SimpleReader) Read(b []byte) (n int, err error) {
	return reader.Read(b)
}

func NewSimpleReader(filename string) *SimpleReader {
	var answer SimpleReader = SimpleReader{
		File: MustOpen(filename),
	}
	switch true {
	case strings.HasSuffix(filename, ".gz"):
		gzipReader, err := gzip.NewReader(answer.File)
		common.ExitIfError(err)
		answer.Reader = bufio.NewReader(gzipReader)
	default:
		answer.Reader = bufio.NewReader(answer.File)
	}
	return &answer
}

func ReadLine(reader *SimpleReader) ([]byte, bool) {
	if curr, err := reader.ReadBytes('\n'); err == nil {
		if curr[len(curr)-1] == '\n' {
			return curr[:len(curr)-1], false
		} else {
			log.Fatalf("Error: end of line did not end with a white space character...\n")
		}
	} else {
		CatchErrThrowEOF(err)
		reader.Close()
	}
	return nil, true
}

func CatchErrThrowEOF(err error) {
	if err == io.EOF {
		return
	} else {
		common.ExitIfError(err)
	}
}

func (reader *SimpleReader) Close() {
	if reader != nil {
		err := reader.File.Close()
		common.ExitIfError(err)
	}
}
