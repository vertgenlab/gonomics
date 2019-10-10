package fastq

import (

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
	file := fileio.EasyOpen(filename)
	defer file.Close()

	answer := ReadFastqs(file)
	return answer
}

func processFastqRecord(line1 string, line2 string, line3 string, line4 string) *Fastq {
	var curr Fastq
	if line3 != "+" {
		log.Fatalf("Error: This line should be a + (plus) sign \n")
	}
	curr = Fastq{Name: line1[1:len(line1)], Seq: dna.StringToBases(line2), Qual: []rune(line4)}
	return &curr
}

func NextFastq(reader *fileio.EasyReader) (*Fastq, bool) {
	line, done := fileio.EasyNextLine(reader)
	line2, done2 := fileio.EasyNextLine(reader)
	line3, done3 := fileio.EasyNextLine(reader)
	line4, done4 := fileio.EasyNextLine(reader)
	if done {
		return nil, true
	}
	if done2 || done3 || done4 {
		log.Fatalf("Error: There is an empty line in this fastq record\n")
	}
	return processFastqRecord(line, line2, line3, line4), false
}

func ReadFastqs(er *fileio.EasyReader) []*Fastq {
	var curr *Fastq
	var done bool
	var answer []*Fastq
	for curr, done = NextFastq(er); !done; curr, done = NextFastq(er) {
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
