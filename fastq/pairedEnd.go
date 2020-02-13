package fastq

import (
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

func PairEndToChan(readOne string, readTwo string, output chan<- *PairedEnd) {
	var curr *PairedEnd
	var done bool

	fileOne := fileio.EasyOpen(readOne)
	defer fileOne.Close()
	fileTwo := fileio.EasyOpen(readTwo)
	defer fileTwo.Close()

	for curr, done = NextFastqPair(fileOne, fileTwo); !done; curr, done = NextFastqPair(fileOne, fileTwo) {
		output <- curr
	}
	close(output)
}

func processFastqPair(readOne *Fastq, readTwo *Fastq) *PairedEnd {
	var curr PairedEnd

	name1 := strings.Split(readOne.Name, " ")
	name2 := strings.Split(readTwo.Name, " ")

	if name1[0] != name2[0] || (string(name1[1][:1]) != "1" && string(name2[1][:1]) != "2") {
		log.Fatalf("Error: Names of fastq pair are not the same\n")
	}
	curr = PairedEnd{Fwd: readOne, Rev: readTwo}
	return &curr
}

func NextFastqPair(reader1 *fileio.EasyReader, reader2 *fileio.EasyReader) (*PairedEnd, bool) {
	fqOne, done1 := NextFastq(reader1)
	fqTwo, done2 := NextFastq(reader2)
	if done1 || done2 {
		return processFastqPair(fqOne, fqTwo), true
	}
	return processFastqPair(fqOne, fqTwo), false
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
