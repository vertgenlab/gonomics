package fastq

import (
	"bufio"
	"fmt"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/fileio"
	//"github.com/vertgenlab/gonomics/qDna"
	"log"
	"math"
)

type Fastq struct {
	Name string
	Seq  []dna.Base
	Qual []rune
}

func Read(filename string) []*Fastq {
	var answer []*Fastq
	var curr *Fastq
	var line string
	var doneName bool = false
	var doneSeq, donePlus, donePhred bool
	var lineSeq, plus, tmpPhred string
	file := fileio.MustOpen(filename)
	defer file.Close()
	reader := bufio.NewReader(file)

	//fastq fields
	var sName string
	var sequence []dna.Base
	var qPhred []rune

	for line, doneName = fileio.NextLine(reader); !doneName; line, doneName = fileio.NextLine(reader) {
		lineSeq, doneSeq = fileio.NextLine(reader)
		plus, donePlus = fileio.NextLine(reader)
		tmpPhred, donePhred = fileio.NextLine(reader)

		if doneSeq || donePlus || donePhred {
			log.Fatalf("Error: lines in %s, must be a multiple of four\n", filename)
		}
		if plus != "+" {
			log.Fatalf("Error: every fourth line in %s should be blank\n", filename)
		}
		sName = line[1:len(line)]
		sequence = dna.StringToBases(lineSeq)
		qPhred = []rune(tmpPhred)
		curr = &Fastq{Name: sName, Seq: sequence, Qual: qPhred}
		answer = append(answer, curr)
	}
	return answer
}

func PhredToPError(ascii rune) float32 {
	q := float64(ascii) - 33
	p := math.Pow(10, -q/10)
	return float32(p)
}

func ErrorRate(ASCII []rune) []float32 {
	var answer []float32
	for i := 0; i < len(ASCII); i++ {
		answer = append(answer, PhredToPError(ASCII[i]))
	}
	return answer
}

func PrintFastq(fq []*Fastq) {
	for i := 0; i < len(fq); i++ {
		readName := "@" + fq[i].Name
		read := dna.BasesToString(fq[i].Seq)
		quality := string(fq[i].Qual)
		fmt.Printf("%s\n%s\n%s\n%s\n", readName, read, "+", quality)
	}
}
