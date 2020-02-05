package fastq

import (
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/fileio"
	"log"
	"strings"
)

type PairedEnd struct {
	Fwd *Fastq
	Rev *Fastq
}

func ReadPairs(readOne string, readTwo string) []*PairedEnd {
	file1 := fileio.EasyOpen(readOne)
	defer file1.Close()

	file2 := fileio.EasyOpen(readTwo)
	defer file2.Close()

	answer := ReadFastqsPairs(file1, file2)
	return answer
}

func processFastqPair(line1 string, line2 string, line3 string, line4 string, line1b string, line2b string, line3b string, line4b string) *PairedEnd {
	var curr PairedEnd

	if line3 != "+" || line3b != "+" {
		log.Fatalf("Error: This line should be a + (plus) sign \n")
	}
	//name := strings.Split(line1, " ")
	//curr = Fastq{Name: name[0][1:len(name[0])], Seq: dna.StringToBases(line2), Qual: []rune(line4)}
	curr = PairedEnd{Fwd: &Fastq{Name: line1[1:len(line1)], Seq: dna.StringToBases(line2), Qual: []rune(line4)}, Rev: &Fastq{Name: line1b[1:len(line1b)], Seq: dna.StringToBases(line2b), Qual: []rune(line4b)}}
	return &curr
}

func NextFastqPair(reader1 *fileio.EasyReader, reader2 *fileio.EasyReader) (*PairedEnd, bool) {
	line1, done := fileio.EasyNextLine(reader1)
	line2, done2 := fileio.EasyNextLine(reader1)
	line3, done3 := fileio.EasyNextLine(reader1)
	line4, done4 := fileio.EasyNextLine(reader1)

	line1b, doneb := fileio.EasyNextLine(reader2)
	line2b, done2b := fileio.EasyNextLine(reader2)
	line3b, done3b := fileio.EasyNextLine(reader2)
	line4b, done4b := fileio.EasyNextLine(reader2)

	if done || doneb {
		return nil, true
	}
	if done2 || done3 || done4 || done2b || done3b || done4b {
		log.Fatalf("Error: There is an empty line in this fastq record\n")
	}

	name1 := strings.Split(line1, " ")
	name2 := strings.Split(line1b, " ")

	if name1[0] != name2[0] || (string(name1[1][:1]) != "1" && string(name2[1][:1]) != "2") {
		log.Fatalf("Error: Names of fastq pair are not the same\n")
	}
	return processFastqPair(line1, line2, line3, line4, line1b, line2b, line3b, line4b), false
}

func ReadFastqsPairs(er *fileio.EasyReader, er2 *fileio.EasyReader) []*PairedEnd {
	var curr *PairedEnd
	var done bool
	var answer []*PairedEnd
	for curr, done = NextFastqPair(er, er2); !done; curr, done = NextFastqPair(er, er2) {
		answer = append(answer, curr)
	}
	return answer
}