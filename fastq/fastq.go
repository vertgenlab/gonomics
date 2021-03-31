package fastq

import (
	"fmt"
	"github.com/vertgenlab/gonomics/common"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/fileio"
	"log"
	"sync"
)

type Fastq struct {
	Name string
	Seq  []dna.Base
	Qual []uint8
}

func Read(filename string) []Fastq {
	file := fileio.EasyOpen(filename)
	var curr Fastq
        var done bool
        var answer []Fastq
        for curr, done = NextFastq(file); !done; curr, done = NextFastq(file) {
                answer = append(answer, curr)
        }
	err := file.Close()
	exception.PanicOnErr(err)
	return answer
}

func ReadToChan(filename string, data chan<- Fastq) {
	file := fileio.EasyOpen(filename)
	for curr, done := NextFastq(file); !done; curr, done = NextFastq(file) {
		data <- curr
	}
	close(data)
	err := file.Close()
	exception.PanicOnErr(err)
}

func GoReadToChan(filename string) <-chan Fastq {
	data := make(chan *Fastq, 1000)
	go ReadToChan(filename, data)
	return data
}

func Write(filename string, records []Fastq) {
	file := fileio.EasyCreate(filename)
	for i := range records {
		WriteToFileHandle(file, records[i])
	}
	err := file.Close()
	exception.PanicOnErr(err)
}

func WriteToFileHandle(file *fileio.EasyWriter, fq Fastq) {
	var err error
	_, err = fmt.Fprintf(file, "%s", ToString(fq))
	exception.PanicOnErr(err)
}

func processFastqRecord(line1 string, line2 string, line3 string, line4 string) Fastq {
	if line3 != "+" {
		log.Fatalf("Error: This line should be a + (plus) sign:\n%s\n", line3)
	}
	curr := Fastq{Name: line1[1:len(line1)], Seq: dna.StringToBases(line2), Qual: ToQual([]byte(line4))}
	return curr
}

func NextFastq(reader *fileio.EasyReader) (Fastq, bool) {
	line, done := fileio.EasyNextLine(reader)
	line2, done2 := fileio.EasyNextLine(reader)
	line3, done3 := fileio.EasyNextLine(reader)
	line4, done4 := fileio.EasyNextLine(reader)
	if done {
		return nil, true
	}
	if done2 || done3 || done4 {
		log.Panicf("Error: Lines in fastq file should be a multiple of 4.\n", )
	}
	return processFastqRecord(line, line2, line3, line4), false
}

func Copy(a *Fastq) *Fastq {
	var answer Fastq = Fastq{Name: a.Name, Seq: make([]dna.Base, len(a.Seq)), Qual: make([]uint8, len(a.Qual))}
	copy(answer.Seq, a.Seq)
	copy(answer.Qual, a.Qual)
	return &answer
}

func ToString(fq Fastq) string {
	readName := "@" + fq[i].Name
	read := dna.BasesToString(fq[i].Seq)
	quality := QualString(fq[i].Qual)
	return fmt.Sprintf("%s\n%s\n%s\n%s\n", readName, read, "+", quality)
}
