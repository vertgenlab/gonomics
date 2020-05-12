package fastq

import (
	"fmt"
	"github.com/vertgenlab/gonomics/common"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/fileio"
	"io"
	"log"
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

func ReadToChan(filename string, output chan<- *Fastq) {
	var curr *Fastq
	var done bool

	file := fileio.EasyOpen(filename)
	defer file.Close()

	for curr, done = NextFastq(file); !done; curr, done = NextFastq(file) {
		output <- curr
	}
	close(output)
}

func Write(filename string, records []*Fastq) {
	file := fileio.EasyCreate(filename)
	defer file.Close()
	WriteToFileHandle(file, records)
}

func WriteToFileHandle(file io.Writer, fq []*Fastq) error {
	var err error
	for i := 0; i < len(fq); i++ {
		_, err = fmt.Fprintf(file, "%s\n%s\n%s\n%s\n", "@"+fq[i].Name, dna.BasesToString(fq[i].Seq), "+", Uint8QualToString(fq[i].Qual))
		common.ExitIfError(err)
	}
	return err
}

func processFastqRecord(line1 string, line2 string, line3 string, line4 string) *Fastq {
	var curr Fastq
	if line3 != "+" {
		log.Fatalf("Error: This line should be a + (plus) sign \n")
	}
	curr = Fastq{Name: line1[1:len(line1)], Seq: dna.StringToBases(line2), Qual: ToQualUint8([]rune(line4))}
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
		quality := string(fq[i].Qual)
		fmt.Printf("%s\n%s\n%s\n%s\n", readName, read, "+", quality)
	}
}
