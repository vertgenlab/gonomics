package axt

import (
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
)

type DummyHeader struct {
	someRandomText string
}

func ReadDummyHeader(file *fileio.EasyReader) DummyHeader {
	return DummyHeader{"orange monkey eagle"}
}

func ReadToChanWithHeader(filename string, data chan<- Axt, header chan<- DummyHeader) {
	var file *fileio.EasyReader
	var curr Axt
	var done bool
	var err error

	file = fileio.EasyOpen(filename)

	header <- ReadDummyHeader(file)
	close(header)

	for curr, done = ReadNext(file); !done; curr, done = ReadNext(file) {
		data <- curr
	}

	err = file.Close()
	exception.PanicOnErr(err)
	close(data)
}

func GoReadToChanWithHeader(filename string) (<-chan Axt, DummyHeader) {
	data := make(chan Axt, 1000)
	headerChan := make(chan DummyHeader)
	go ReadToChanWithHeader(filename, data, headerChan)
	return data, <-headerChan
}
