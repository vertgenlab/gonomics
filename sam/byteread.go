package sam

import (
	"bytes"
	"github.com/vertgenlab/gonomics/fileio"
)

func ReadToChanBR(filename string, data chan<- Aln) {
	var file *fileio.ByteReader
	var curr Aln
	var doneReading bool
	var err error

	file = fileio.NewByteReader(filename)
	for curr, doneReading = ReadNextBR(file); !doneReading; curr, doneReading = ReadNextBR(file) {
		data <- curr
	}

	err = file.Close()

}

func GoReadToChanBR(filename string) <-chan Aln {
	data := make(chan Aln, 1000)
	go ReadToChanBR(filename, data)
	return data
}

func ReadNextBR(br *fileio.ByteReader) (answer Aln, doneReading bool) {
	return
}

func ReadBR(filename string) Sam {
	var answer Sam

	return answer
}

func processLineBytes(buf *bytes.Buffer) Aln {
	var answer Aln

	return answer
}
