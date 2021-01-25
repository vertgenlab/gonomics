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

func Read(filename string) []*Fastq {
	file := fileio.EasyOpen(filename)
	defer file.Close()

	answer := ReadFastqs(file)
	return answer
}

func ReadToChan(file *fileio.EasyReader, data chan<- *Fastq, wg *sync.WaitGroup) {
	for curr, done := NextFastq(file); !done; curr, done = NextFastq(file) {
		data <- curr
	}
	file.Close()
	wg.Done()
}

func GoReadToChan(filename string) <-chan *Fastq {
	file := fileio.EasyOpen(filename)
	var wg sync.WaitGroup
	data := make(chan *Fastq)
	wg.Add(1)
	go ReadToChan(file, data, &wg)

	go func() {
		wg.Wait()
		close(data)
	}()

	return data
}

func Write(filename string, records []*Fastq) {
	file := fileio.EasyCreate(filename)
	defer file.Close()
	for _, fq := range records {
		WriteToFileHandle(file, fq)
	}
}

func WriteToFileHandle(file *fileio.EasyWriter, fq *Fastq) {
	var err error
	_, err = fmt.Fprintf(file, "@%s\n%s\n+\n%s\n", fq.Name, dna.BasesToString(fq.Seq), QualString(fq.Qual))
	common.ExitIfError(err)
}

func processFastqRecord(line1 string, line2 string, line3 string, line4 string) *Fastq {
	var curr *Fastq
	if line3 != "+" {
		log.Fatalf("Error: This line should be a + (plus) sign \n")
	}
	currBases, err := dna.StringToBases(line2)
	if err != nil {
		log.Panicf("error converting to bases")
	}
	curr = &Fastq{Name: line1[1:len(line1)], Seq: currBases, Qual: ToQual([]byte(line4))}
	return curr
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

func Copy(a *Fastq) *Fastq {
	var answer Fastq = Fastq{Name: a.Name, Seq: make([]dna.Base, len(a.Seq)), Qual: make([]uint8, len(a.Qual))}
	copy(answer.Seq, a.Seq)
	copy(answer.Qual, a.Qual)
	return &answer
}

func PrintFastq(fq []*Fastq) {
	for i := 0; i < len(fq); i++ {
		readName := "@" + fq[i].Name
		read := dna.BasesToString(fq[i].Seq)
		quality := QualString(fq[i].Qual)
		fmt.Printf("%s\n%s\n%s\n%s\n", readName, read, "+", quality)
	}
}
