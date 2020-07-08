package fasta

import (
	"fmt"
	"github.com/vertgenlab/gonomics/common"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/fileio"
	"io"
	"log"
	"os"
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

func WriteToSplitChr(filename string, records []*Fasta) {
	for _, rec := range records {
		Write(filename+rec.Name+".fa", []*Fasta{rec})
	}
}

func WriteToFileHandle(file io.Writer, records []*Fasta, lineLength int) {
	var err error
	for _, rec := range records {
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
}

func Write(filename string, records []*Fasta) {
	lineLength := 50
	file := fileio.EasyCreate(filename)
	defer file.Close()

	WriteToFileHandle(file, records, lineLength)
}

func WriteFasta(file io.Writer, rec *Fasta, lineLength int) {
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
	//New line character
	_, err = fmt.Fprint(file, "\n")
	common.ExitIfError(err)
}

func WriteGroups(filename string, groups [][]*Fasta) error {
	lineLength := 50
	file, err := os.Create(filename)
	if err != nil {
		return err
	}
	defer file.Close()

	for i, _ := range groups {
		WriteToFileHandle(file, groups[i], lineLength)
		_, err = fmt.Fprint(file, "\n")
		common.ExitIfError(err)
	}
	return nil
}

func CreateAllGaps(name string, numGaps int64) *Fasta {
	answer := Fasta{Name: name, Seq: dna.CreateAllGaps(numGaps)}
	return &answer
}

func CreateAllNs(name string, numGaps int64) *Fasta {
	answer := Fasta{Name: name, Seq: dna.CreateAllNs(numGaps)}
	return &answer
}

//Dictionary/hash map look up of sequence by name
func FastaMap(ref []*Fasta) map[string][]dna.Base {
	m := make(map[string][]dna.Base)
	var curr *Fasta
	for i := 0; i < len(ref); i++ {
		curr = ref[i]
		_, ok := m[curr.Name]
		if !ok {
			m[curr.Name] = curr.Seq
		} else {
			log.Fatalf("Fasta slice has duplicate names")
		}
	}
	return m
}
