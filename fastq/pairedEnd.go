package fastq

import (
	"fmt"
	"github.com/vertgenlab/gonomics/common"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/fileio"
	"log"
	"strings"
	"sync"
)

type PairedEnd struct {
	Fwd *Fastq
	Rev *Fastq
}

type PairedEndBig struct {
	Fwd FastqBig
	Rev FastqBig
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

func ReadPairBigToChan(readOne string, readTwo string, answer chan<- PairedEndBig) {
	var fq PairedEndBig
	fileOne, fileTwo := fileio.EasyOpen(readOne), fileio.EasyOpen(readTwo)

	defer fileOne.Close()
	defer fileTwo.Close()

	for curr, done := NextFastqPair(fileOne, fileTwo); !done; curr, done = NextFastqPair(fileOne, fileTwo) {
		fq = PairedEndBig{Fwd: *ToFastqBig(curr.Fwd), Rev: *ToFastqBig(curr.Rev)}
		answer <- fq
	}
	close(answer)
}

func NextFastqPair(reader1 *fileio.EasyReader, reader2 *fileio.EasyReader) (*PairedEnd, bool) {
	fqOne, done1 := NextFastq(reader1)
	fqTwo, done2 := NextFastq(reader2)
	if (!done1 && done2) || (done1 && !done2) {
		log.Fatalf("Error: fastq files do not end at the same time...\n")
	} else if done1 || done2 {
		return nil, true
	}
	curr := &PairedEnd{Fwd: nil, Rev: nil}
	curr.Fwd, curr.Rev = fqOne, fqTwo
	curr.Fwd.Name, curr.Rev.Name = strings.Split(fqOne.Name, " ")[0], strings.Split(fqTwo.Name, " ")[0]
	return curr, false
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

func WritingHelper(fileOne *fileio.EasyWriter, fileTwo *fileio.EasyWriter, fq *PairedEnd) {
	//TODO: figure out why this seems a little slower
	//WriteToFileHandle(fileOne, fq.Fwd)
	//WriteToFileHandle(fileTwo, fq.Rev)
	var err error
	_, err = fmt.Fprintf(fileOne, "@%s\n%s\n+\n%s\n", fq.Fwd.Name, dna.BasesToString(fq.Fwd.Seq), QualString(fq.Fwd.Qual))
	common.ExitIfError(err)
	_, err = fmt.Fprintf(fileTwo, "@%s\n%s\n+\n%s\n", fq.Rev.Name, dna.BasesToString(fq.Rev.Seq), QualString(fq.Rev.Qual))
	common.ExitIfError(err)
}

func WritePair(readOne string, readTwo string, records []*PairedEnd) {
	fileOne := fileio.EasyCreate(readOne)
	fileTwo := fileio.EasyCreate(readTwo)
	defer fileOne.Close()
	defer fileTwo.Close()
	for _, fq := range records {
		WritingHelper(fileOne, fileTwo, fq)
	}
}

func WritingChan(readOne string, readTwo string, output <-chan *PairedEnd, wg *sync.WaitGroup) {
	fileOne, fileTwo := fileio.EasyCreate(readOne), fileio.EasyCreate(readTwo)
	defer fileOne.Close()
	defer fileTwo.Close()
	for fq := range output {
		WritingHelper(fileOne, fileTwo, fq)
	}
	wg.Done()
}

func GoWriteFqPair(readOne string, readTwo string, data <-chan *PairedEnd) {
	var wg sync.WaitGroup
	wg.Add(1)
	go WritingChan(readOne, readTwo, data, &wg)
	wg.Wait()
}

// ReadFqBigPair will take 2 readers which will convert paired end read fastq files into a paired end FastqBig struct.
// Note: while this function is return as a pointer, it's purpose is to be deferenced at the next function call.
// In additon the pointers to read one and read two will also be dereference. When sending a pointer to a struct through
// a channel, or a struct with pointers inside, memory allocated will be placed on the heap hindering performance.
func ReadFqBigPair(readerOne *fileio.SimpleReader, readerTwo *fileio.SimpleReader) (*PairedEndBig, bool) {
	var done bool
	var fqOne, fqTwo *FastqBig = &FastqBig{}, &FastqBig{}
	fqOne, done = ReadFqBig(readerOne)
	if !done || fqOne != nil {
		fqTwo, done = ReadFqBig(readerTwo)
		if !done || fqTwo != nil {
			return &PairedEndBig{Fwd: *fqOne, Rev: *fqTwo}, false
		} else {
			log.Fatalf("Error: fastq files do not end at the same time...\n")
		}
	}
	return nil, true
}
