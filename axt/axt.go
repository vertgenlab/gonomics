// Package axt provides the struct and functions that operate on alignments in Axt format.
package axt

import (
	"fmt"
	"github.com/vertgenlab/gonomics/common"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"io"
	"log"
	"strings"
)

// Axt struct: Naming convention is hard here because UCSC website does not
// match the UCSC Kent source tree.
type Axt struct {
	RName      string
	RStart     int
	REnd       int
	QName      string
	QStart     int
	QEnd       int
	QStrandPos bool // true is positive strand, false is negative strand
	Score      int
	RSeq       []dna.Base
	QSeq       []dna.Base
}

// Read takes a filename and returns a slice of all Axt records found in the file.
func Read(filename string) []Axt {
	var answer []Axt
	var curr Axt
	var done bool
	var err error
	var file *fileio.EasyReader

	file = fileio.EasyOpen(filename)
	for curr, done = ReadNext(file); !done; curr, done = ReadNext(file) {
		answer = append(answer, curr)
	}
	err = file.Close()
	exception.PanicOnErr(err)
	return answer
}

// ReadToChan takes a filename, a channel for the axt structs, and a channel for
// the header text.  The function will read all remaining Axt records
// from the open file into the channel.  The function will then close the file and close the channel.
func ReadToChan(filename string, data chan<- Axt, header chan<- []string) {
	var file *fileio.EasyReader
	var curr Axt
	var done bool
	var err error
	var headerLines []string

	file = fileio.EasyOpen(filename)

	headerLines, err = fileio.EasyReadHeader(file)
	exception.PanicOnErr(err)
	header <- headerLines
	close(header)

	for curr, done = ReadNext(file); !done; curr, done = ReadNext(file) {
		data <- curr
	}

	err = file.Close()
	exception.PanicOnErr(err)
	close(data)
}

// GoReadToChan takes a filename and returns a channel and header text.  GoReadToChan
// will launch a Go routine to read all the alignments from the file into the channel
// and then close the channel.
func GoReadToChan(filename string) (<-chan Axt, []string) {
	data := make(chan Axt, 1000)
	header := make(chan []string)
	go ReadToChan(filename, data, header)
	return data, <-header
}

// ReadNext takes an EasyReader and returns the next Axt record as well as a boolean flag
// indicating if we are done reading the file.  If there
// is an Axt record left in the file, it will be returned along with "false."  If we are done reading the file
// a blank Axt record will be returned along with "true."
func ReadNext(reader *fileio.EasyReader) (Axt, bool) {
	header, hDone := fileio.EasyNextRealLine(reader)
	if hDone {
		return Axt{}, true
	}
	rSeq, rDone := fileio.EasyNextRealLine(reader)
	qSeq, qDone := fileio.EasyNextRealLine(reader)
	blank, bDone := fileio.EasyNextRealLine(reader)
	if blank != "" {
		log.Fatalf("Error: every fourth line in an axt file should be blank\n")
	}
	if rDone || qDone || bDone {
		log.Fatalf("Error: number of lines in an axt file must be a multiple of four\n")
	}
	return axtHelper(header, rSeq, qSeq), false
}

// axtHelper is a helper function to process individual axt records.
func axtHelper(header string, rSeq string, qSeq string) Axt {
	var words []string = strings.Split(header, " ")
	if len(words) != 9 {
		log.Fatalf("Error: expecting 9 space-separated fields in the first line of each axt record.  Found:%d in line:\n%s\n", len(words), header)
	}
	if rSeq == "" || qSeq == "" {
		log.Fatalf("Error: missing reference seq or query seq in axt record:\n%s\n", header)
	}
	var answer Axt = Axt{
		RName:      words[1],
		RStart:     common.StringToInt(words[2]),
		REnd:       common.StringToInt(words[3]),
		QName:      words[4],
		QStart:     common.StringToInt(words[5]),
		QEnd:       common.StringToInt(words[6]),
		QStrandPos: common.StringToStrand(words[7]),
		Score:      common.StringToInt(words[8]),
		RSeq:       dna.StringToBases(rSeq),
		QSeq:       dna.StringToBases(qSeq),
	}
	return answer
}

// WriteToFileHandle writes a given axt record to a file.  Axt records are numbered in files, so
// a number must also be given.  These numbers usually start at zero and count up throughout the file.
func WriteToFileHandle(file io.Writer, input Axt, alnNumber int) {
	_, err := fmt.Fprintf(file, "%s", ToString(input, alnNumber))
	exception.PanicOnErr(err)
}

// ToString converts an Axt alignment struct into a string.  Axt records are numbered in files, so
// a number must be provided that will be used as this id.
func ToString(input Axt, id int) string {
	return fmt.Sprintf("%d %s %d %d %s %d %d %c %d\n%s\n%s\n\n", id, input.RName, input.RStart, input.REnd, input.QName, input.QStart, input.QEnd, common.StrandToRune(input.QStrandPos), input.Score, dna.BasesToString(input.RSeq), dna.BasesToString(input.QSeq))
}

// Write is a wrapper function that will loop over a slice of axt alignments and writes each record to a file.
func Write(filename string, data []Axt) {
	var file *fileio.EasyWriter
	var err error
	var i int

	file = fileio.EasyCreate(filename)
	for i = range data {
		WriteToFileHandle(file, data[i], i)
	}
	err = file.Close()
	exception.PanicOnErr(err)
}

// Swap will swap reference and query in the alignment.
func Swap(in *Axt, tLen int, qLen int) {
	in.RSeq, in.QSeq = in.QSeq, in.RSeq
	in.RName, in.QName = in.QName, in.RName
	if !in.QStrandPos {
		in.RStart, in.REnd = qLen-in.QEnd+1, qLen-in.QStart+1
		in.QStart, in.QEnd = tLen-in.REnd+1, tLen-in.RStart+1
		dna.ReverseComplement(in.RSeq)
		dna.ReverseComplement(in.QSeq)
	} else {
		in.RStart, in.REnd = in.QStart, in.QEnd
		in.QStart, in.QEnd = in.RStart, in.REnd
	}
	in.RSeq, in.QSeq = in.QSeq, in.RSeq
}

// IsAxtFile returns true of filename ends in .axt or .axt.gz
func IsAxtFile(filename string) bool {
	if strings.HasSuffix(filename, ".axt") || strings.HasSuffix(filename, ".axt.gz") {
		return true
	} else {
		return false
	}
}
