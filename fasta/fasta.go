package fasta

import (
	"fmt"
	"github.com/vertgenlab/gonomics/common"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/fileio"
	"strings"
)

type Fasta struct {
	Name string
	Seq  []dna.Base
}

func Read(filename string) []*Fasta {
	var line string
	var currSeq []dna.Base
	var answer []*Fasta
	var seqIdx int64 = -1
	var doneReading bool = false

	file := fileio.EasyOpen(filename)
	defer file.Close()
	//reader := bufio.NewReader(file)

	for line, doneReading = fileio.EasyNextRealLine(file); !doneReading; line, doneReading = fileio.EasyNextRealLine(file) {
		if strings.HasPrefix(line, ">") {
			tmp := Fasta{Name: line[1:len(line)]}
			answer = append(answer, &tmp)
			seqIdx++
		} else {
			currSeq = dna.StringToBases(line)
			answer[seqIdx].Seq = append(answer[seqIdx].Seq, currSeq...)
		}
	}
	return answer
}

func WriteToFileHandle(file *fileio.EasyWriter, rec *Fasta, lineLength int) {
	var err error
	_, err = fmt.Fprintf(file, ">%s\n", rec.Name)
	common.ExitIfError(err)
	for i := 0; i < len(rec.Seq); i += lineLength {
		if i+lineLength > len(rec.Seq) {
			_, err = fmt.Fprintf(file, "%s\n", dna.BasesToString(rec.Seq[i:]))
			common.ExitIfError(err)
		} else {
			_, err = fmt.Fprintf(file, "%s\n", dna.BasesToString(rec.Seq[i:i+lineLength]))
			common.ExitIfError(err)
		}
	}
}

func Write(filename string, records []*Fasta) {
	lineLength := 50
	file := fileio.EasyCreate(filename)
	defer file.Close()
	for _, rec := range records {
		WriteToFileHandle(file, rec, lineLength)
	}
}

func CreateAllGaps(name string, numGaps int64) *Fasta {
	answer := Fasta{Name: name, Seq: dna.CreateAllGaps(numGaps)}
	return &answer
}

func CreateAllNs(name string, numGaps int64) *Fasta {
	answer := Fasta{Name: name, Seq: dna.CreateAllNs(numGaps)}
	return &answer
}
