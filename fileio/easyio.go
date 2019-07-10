package fileio

import (
	"bufio"
	"compress/gzip"
	"github.com/vertgenlab/gonomics/common"
	"os"
	"strings"
)

type EasyReader struct {
	file         *os.File
	internalBuff *bufio.Reader
	internalGzip *gzip.Reader
	BuffReader   *bufio.Reader
}

type EasyWriter struct {
	file         *os.File
	internalBuff *bufio.Writer
	internalGzip *gzip.Writer
}

func EasyOpen(filename string) *EasyReader {
	answer := EasyReader{}
	answer.file = MustOpen(filename)
	var err error

	if strings.HasSuffix(filename, ".gz") {
		answer.internalBuff = bufio.NewReader(answer.file)
		answer.internalGzip, err = gzip.NewReader(answer.internalBuff)
		common.ExitIfError(err)
		answer.BuffReader = bufio.NewReader(answer.internalGzip)
	} else {
		answer.BuffReader = bufio.NewReader(answer.file)
		answer.internalBuff = nil
		answer.internalGzip = nil
	}
	return &answer
}

func EasyCreate(filename string) *EasyWriter {
	answer := EasyWriter{}
	answer.file = MustCreate(filename)
	answer.internalBuff = bufio.NewWriter(answer.file)

	if strings.HasSuffix(filename, ".gz") {
		answer.internalGzip = gzip.NewWriter(answer.internalBuff)
	} else {
		answer.internalGzip = nil
	}
	return &answer
}

func EasyNextLine(file *EasyReader) (string, bool) {
	return NextLine(file.BuffReader)
}

func EasyNextRealLine(file *EasyReader) (string, bool) {
	return NextRealLine(file.BuffReader)
}

func (er *EasyReader) Close() {
	/*if er.BuffReader != nil {
		er.BuffReader.Close()
	}*/
	if er.internalGzip != nil {
		er.internalGzip.Close()
	}
	/*if er.internalBuff != nil {
		er.internalBuff.Close()
	}*/
	if er.file != nil {
		er.file.Close()
	}
}

func (ew *EasyWriter) Close() {
	if ew.internalGzip != nil {
		ew.internalGzip.Close()
	}
	if ew.internalBuff != nil {
		ew.internalBuff.Flush()
	}
	if ew.file != nil {
		ew.file.Close()
	}
}

func (ew *EasyWriter) Write(p []byte) (n int, err error) {
	if ew.internalGzip != nil {
		return ew.internalGzip.Write(p)
	} else {
		return ew.internalBuff.Write(p)
	}
}
