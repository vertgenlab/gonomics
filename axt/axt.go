package axt

import (
	"bufio"
	"fmt"
	"github.com/vertgenlab/gonomics/common"
	"github.com/vertgenlab/gonomics/dna"
	"io"
	"os"
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
	var err, hErr, rErr, qErr, bErr, startErr, endErr error
	var words []string

	file := common.MustOpen(filename)
	defer file.Close()
	reader := bufio.NewReader(file)

	for header, hErr = common.NextRealLine(reader); hErr == nil; header, hErr = common.NextRealLine(reader) {
		rSeq, rErr = common.NextRealLine(reader)
		qSeq, qErr = common.NextRealLine(reader)
		blank, bErr = common.NextRealLine(reader)
		if rErr == io.EOF || qErr == io.EOF || bErr == io.EOF {
			common.Exit(fmt.Sprintf("Error: lines in %s, must be a multiple of four\n", filename))
		}
		if blank != "" {
			common.Exit(fmt.Sprintf("Error: every fourth line in %s should be blank\n", filename))
		}

		words = strings.Split(header, " ")
		if len(words) != 9 {
			common.Exit(fmt.Sprintf("Error: sequences in %s should be the same length\n", header))
		}

		curr := Axt{}
		curr.RName = words[1]
		curr.RStart, startErr = strconv.ParseInt(words[2], 10, 64)
		curr.REnd, endErr = strconv.ParseInt(words[3], 10, 64)
		if startErr != nil || endErr != nil {
			common.Exit(fmt.Sprintf("Error: trouble parsing reference start and end in %s\n", header))
		}
		curr.QName = words[4]
		curr.QStart, startErr = strconv.ParseInt(words[5], 10, 64)
		curr.QEnd, endErr = strconv.ParseInt(words[6], 10, 64)
		if startErr != nil || endErr != nil {
			common.Exit(fmt.Sprintf("Error: trouble parsing query start and end in %s\n", header))
		}
		switch words[7] {
		case "+":
			curr.QStrandPos = true
		case "-":
			curr.QStrandPos = false
		default:
			common.Exit(fmt.Sprintf("Error: did not recognize strand in %s\n", header))
		}
		curr.Score, err = strconv.ParseInt(words[8], 10, 64)
		if err != nil {
			common.Exit(fmt.Sprintf("Error: trouble parsing the score in %s\n", header))
		}
		curr.RSeq, err = dna.StringToBases(rSeq)
		common.ExitIfError(err)
		curr.QSeq, err = dna.StringToBases(qSeq)
		common.ExitIfError(err)

		answer = append(answer, &curr)
	}
	if hErr != io.EOF {
		common.ExitIfError(hErr)
	}
	return answer
}

func WriteToFileHandle(file *os.File, input *Axt, alnNumber int) error {
	var strand rune
	if input.QStrandPos {
		strand = '+'
	} else {
		strand = '-'
	}
	_, err := fmt.Fprintf(file, "%d %s %d %d %s %d %d %c %d\n%s\n%s\n\n", alnNumber, input.RName, input.RStart, input.REnd, input.QName, input.QStart, input.QEnd, strand, input.Score, dna.BasesToString(input.RSeq), dna.BasesToString(input.QSeq))
	return err
}

func Write(filename string, data []*Axt) error {
	file, err := os.Create(filename)
	if err != nil {
		return err
	}
	defer file.Close()

	for i, _ := range data {
		err = WriteToFileHandle(file, data[i], i)
		if err != nil {
			return err
		}
	}
	return err
}
