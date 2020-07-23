//Package axt provides the struct and functions that operate on axt alignment formats.
package axt

import (
	"fmt"
	"github.com/vertgenlab/gonomics/common"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/fileio"
	"log"
	"strconv"
	"strings"
)

// Axt struct: Naming convention is hard here because UCSC website does not
// match the UCSC Kent source tree.
type Axt struct {
	RName      string
	RStart     int64
	REnd       int64
	QName      string
	QStart     int64
	QEnd       int64
	QStrandPos bool // true is positive strand, false is negative strand
	Score      int64
	RSeq       []dna.Base
	QSeq       []dna.Base
}

//Read is a function enabling the processing of axt text files into a data structure used by gonomics.
func Read(filename string) []*Axt {
	var answer []*Axt
	var header, rSeq, qSeq, blank string
	var err, startErr, endErr error
	var hDone, rDone, qDone, bDone bool
	var words []string

	file := fileio.EasyOpen(filename)
	defer file.Close()
	for header, hDone = fileio.EasyNextRealLine(file); !hDone; header, hDone = fileio.EasyNextRealLine(file) {
		rSeq, rDone = fileio.EasyNextRealLine(file)
		qSeq, qDone = fileio.EasyNextRealLine(file)
		blank, bDone = fileio.EasyNextRealLine(file)
		if rDone || qDone || bDone {
			log.Fatalf("Error: lines in %s, must be a multiple of four\n", filename)
		}
		if blank != "" {
			log.Fatalf("Error: every fourth line in %s should be blank\n", filename)
		}

		words = strings.Split(header, " ")
		if len(words) != 9 {
			log.Fatalf("Error: sequences in %s should be the same length\n", header)
		}

		curr := Axt{}
		curr.RName = words[1]
		curr.RStart, startErr = strconv.ParseInt(words[2], 10, 64)
		curr.REnd, endErr = strconv.ParseInt(words[3], 10, 64)
		if startErr != nil || endErr != nil {
			log.Fatalf("Error: trouble parsing reference start and end in %s\n", header)
		}
		curr.QName = words[4]
		curr.QStart, startErr = strconv.ParseInt(words[5], 10, 64)
		curr.QEnd, endErr = strconv.ParseInt(words[6], 10, 64)
		if startErr != nil || endErr != nil {
			log.Fatalf("Error: trouble parsing query start and end in %s\n", header)
		}
		switch words[7] {
		case "+":
			curr.QStrandPos = true
		case "-":
			curr.QStrandPos = false
		default:
			log.Fatalf("Error: did not recognize strand in %s\n", header)
		}
		curr.Score, err = strconv.ParseInt(words[8], 10, 64)
		if err != nil {
			log.Fatalf("Error: trouble parsing the score in %s\n", header)
		}
		curr.RSeq = dna.StringToBases(rSeq)
		curr.QSeq = dna.StringToBases(qSeq)

		answer = append(answer, &curr)
	}
	return answer
}

//ReadToChan is a function that takes an EasyReader, which uses type os.File as an input to read axt alignments into an Axt channel.
func ReadToChan(reader *fileio.EasyReader, answer chan<- *Axt) {
	for data, err := NextAxt(reader); !err; data, err = NextAxt(reader) {
		answer <- data
	}
	close(answer)
}

func GoReadToChan(filename string) chan *Axt {
	answer := make(chan *Axt)
	file := fileio.EasyOpen(filename)
	go ReadToChan(file, answer)
	return answer
}

//NextAxt processes the next Axt alignment in the provided input.
func NextAxt(reader *fileio.EasyReader) (*Axt, bool) {
	header, hDone := fileio.EasyNextRealLine(reader)
	rSeq, rDone := fileio.EasyNextRealLine(reader)
	qSeq, qDone := fileio.EasyNextRealLine(reader)
	blank, bDone := fileio.EasyNextRealLine(reader)
	if blank != "" {
		log.Fatalf("Error: every fourth line should be blank\n")
	}
	if hDone || rDone || qDone || bDone {
		return nil, true
	}
	return axtHelper(header, rSeq, qSeq, blank), false
}

//axtHelper is a helper function to process individual axt alignments.
func axtHelper(header string, rSeq string, qSeq string, blank string) *Axt {
	var words []string = strings.Split(header, " ")
	if len(words) != 9 || rSeq == "" || qSeq == "" {
		log.Fatalf("Error: missing fields in header or sequences\n")
	}
	var answer *Axt = &Axt{
		RName:      words[1],
		RStart:     common.StringToInt64(words[2]),
		REnd:       common.StringToInt64(words[3]),
		QName:      words[4],
		QStart:     common.StringToInt64(words[5]),
		QEnd:       common.StringToInt64(words[6]),
		QStrandPos: common.StringToStrand(words[7]),
		Score:      common.StringToInt64(words[8]),
		RSeq:       dna.StringToBases(rSeq),
		QSeq:       dna.StringToBases(qSeq),
	}
	return answer
}

//WriteToFileHandle writes a given axt record to file as well as handling possible errors.
func WriteToFileHandle(file *fileio.EasyWriter, input *Axt, alnNumber int) {
	_, err := fmt.Fprintf(file, "%s", ToString(input, alnNumber))
	common.ExitIfError(err)
}

//ToString converts an Axt alignment struct into a string.
func ToString(input *Axt, id int) string {
	return fmt.Sprintf("%d %s %d %d %s %d %d %c %d\n%s\n%s\n\n", id, input.RName, input.RStart, input.REnd, input.QName, input.QStart, input.QEnd, common.StrandToRune(input.QStrandPos), input.Score, dna.BasesToString(input.RSeq), dna.BasesToString(input.QSeq))
}

//Write is a wrapper function that will loop over a slice of axt alignments and writes each record to a file.
func Write(filename string, data []*Axt) {
	file := fileio.EasyCreate(filename)
	defer file.Close()

	for i, _ := range data {
		WriteToFileHandle(file, data[i], i)
	}
}

//AxtInfo is a pretty print function that will return a string that contains only the header info from axt record without sequences from target and/or query.
func AxtInfo(input *Axt) string {
	var text string = ""
	text = fmt.Sprintf("%s;%d;%d;%s;%d;%d;%t;%d", input.RName, input.RStart, input.REnd, input.QName, input.QStart, input.QEnd, input.QStrandPos, input.Score)
	return text
}

//SwapBoth will preform a simple swap with target and query records contained inside axt alignment.
func SwapBoth(in *Axt, tLen int64, qLen int64) *Axt {
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
	return in
}

//IsAxtFile checks suffix of file name to confirm axt format
func IsAxtFile(filename string) bool {
	if strings.HasSuffix(filename, ".axt") || strings.HasSuffix(filename, ".axt.gz") {
		return true
	} else {
		return false
	}
}
