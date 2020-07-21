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

// Naming convention is hard here because UCSC website does not
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

func WriteToFileHandle(file *fileio.EasyWriter, input *Axt, alnNumber int) {
	_, err := fmt.Fprintf(file, "%s", ToString(input, alnNumber))
	common.ExitIfError(err)
}

func ToString(input *Axt, id int) string {
	return fmt.Sprintf("%d %s %d %d %s %d %d %c %d\n%s\n%s\n\n", id, input.RName, input.RStart, input.REnd, input.QName, input.QStart, input.QEnd, common.StrandToRune(input.QStrandPos), input.Score, dna.BasesToString(input.RSeq), dna.BasesToString(input.QSeq))
}

func Write(filename string, data []*Axt) {
	file := fileio.EasyCreate(filename)
	defer file.Close()

	for i, _ := range data {
		WriteToFileHandle(file, data[i], i)
	}
}

func AxtInfo(input *Axt) string {
	var text string = ""
	text = fmt.Sprintf("%s;%d;%d;%s;%d;%d;%t;%d", input.RName, input.RStart, input.REnd, input.QName, input.QStart, input.QEnd, input.QStrandPos, input.Score)
	return text
}

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
