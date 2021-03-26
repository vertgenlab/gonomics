package sam

import (
	"bytes"
	"github.com/vertgenlab/gonomics/exception"
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
	exception.PanicOnErr(err)
	close(data)
}

// GoReadToChanBR streams the input file so that only a small portion
// of the file is kept in memory at a time.
func GoReadToChanBR(filename string) <-chan Aln {
	data := make(chan Aln, 1000)
	go ReadToChanBR(filename, data)
	val := <-data
	return data
}

// ReadBR reads the entire file into a Sam struct where each record
// is an index in Sam.Aln. Note that sam files can get very large
// such that storing the entire file in memory is not feasible. Most
// sam files should be read using the GoReadToChan function which
// streams sam records so only a small portion of the file is kept
// in memory at any given time.
func ReadBR(filename string) Sam {
	var answer Sam
	var file *fileio.ByteReader
	var curr Aln
	var doneReading bool
	var err error
	file = fileio.NewByteReader(filename)
	answer.Header = ReadHeaderBR(file)
	for curr, doneReading = ReadNextBR(file); !doneReading; curr, doneReading = ReadNextBR(file) {
		answer.Aln = append(answer.Aln, curr)
	}
	err = file.Close()
	exception.PanicOnErr(err)
	return answer
}

// ReadNext takes a ByteReader and returns the next Sam record as well as a boolean flag
// indicating if the file is finished being read. If there is a Sam record to process
// the function will return the Sam record and 'false'. After processing all Sam records
// in the file, the function will return a blank Aln and 'true'.
func ReadNextBR(br *fileio.ByteReader) (answer Aln, doneReading bool) {
	return
}

// processLineBytes parses a bytes.Buffer filled with a single line from an input Sam file.
func processLineBytes(buf *bytes.Buffer) Aln {
	var answer Aln

	return answer
}
