package fastq

import (
	"fmt"
	"log"
	"strings"
	"sync"

	"github.com/vertgenlab/gonomics/common"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/fileio"
)

// PairedEnd is a struct that contains two paired Fastq structs, marked Fwd and Rev.
type PairedEnd struct {
	Fwd Fastq
	Rev Fastq
}

// PairedEndBig is a struct that contains two paired FastqBig structs, marked Fwd and Rev.
type PairedEndBig struct {
	Fwd FastqBig
	Rev FastqBig
}

// ReadPairs takes two input file names and returns a slice of PairedEnd structs. For large input files, it is advisable to use PairEndToChan.
func ReadPairs(readOne string, readTwo string) []PairedEnd {
	file1 := fileio.EasyOpen(readOne)
	defer file1.Close()

	file2 := fileio.EasyOpen(readTwo)
	defer file2.Close()

	answer := ReadFastqPairs(file1, file2)
	return answer
}

// PairEndToChan parses PairedEnd structs from two input fastq files, with filenames readOne and readTwo, and sends them to a channel named output.
func PairedEndToChan(readOne string, readTwo string, output chan<- PairedEnd) {
	var curr PairedEnd
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

// ReadPairBigToChan is an analog of PairEndToChan for producing input channels of PairedEndBig structs.
func ReadPairBigToChan(fileOne string, fileTwo string, answer chan<- PairedEndBig) {
	readOne, readTwo := fileio.NewByteReader(fileOne), fileio.NewByteReader(fileTwo)
	for fq, done := ReadFqBigPair(readOne, readTwo); !done; fq, done = ReadFqBigPair(readOne, readTwo) {
		answer <- fq
	}
	close(answer)
}

// NextFastqPair is a helper funcotion of PairEndToChan. It checks a reader for additional data lines and parses PairedEnd structs if more lines exist. The bool return indicates 'done', or that a file has no additional data lines.
func NextFastqPair(reader1 *fileio.EasyReader, reader2 *fileio.EasyReader) (PairedEnd, bool) {
	fqOne, done1 := NextFastq(reader1)
	fqTwo, done2 := NextFastq(reader2)
	if (!done1 && done2) || (done1 && !done2) {
		log.Fatalf("Error: fastq files do not end at the same time...\n")
	} else if done1 || done2 {
		return PairedEnd{}, true
	}
	curr := PairedEnd{Fwd: Fastq{}, Rev: Fastq{}}
	curr.Fwd, curr.Rev = fqOne, fqTwo
	curr.Fwd.Name, curr.Rev.Name = strings.Split(fqOne.Name, " ")[0], strings.Split(fqTwo.Name, " ")[0]
	return curr, false
}

// ReadFastqPairs accepts two fileio.EasyReaders and parses PairedEnd structs, returning these structs in a slice of pointers.
func ReadFastqPairs(er *fileio.EasyReader, er2 *fileio.EasyReader) []PairedEnd {
	var curr PairedEnd
	var done bool
	var answer []PairedEnd
	for curr, done = NextFastqPair(er, er2); !done; curr, done = NextFastqPair(er, er2) {
		answer = append(answer, curr)
	}
	return answer
}

// WritingHelper is a helper function of PairedEnd write functions, and converts PairedEnd structs to strings and sends them to the respective writers.
func WritingHelper(fileOne *fileio.EasyWriter, fileTwo *fileio.EasyWriter, fq PairedEnd) {
	//TODO: figure out why this seems a little slower
	//WriteToFileHandle(fileOne, fq.Fwd)
	//WriteToFileHandle(fileTwo, fq.Rev)
	var err error
	_, err = fmt.Fprintf(fileOne, "@%s\n%s\n+\n%s\n", fq.Fwd.Name, dna.BasesToString(fq.Fwd.Seq), QualString(fq.Fwd.Qual))
	common.ExitIfError(err)
	_, err = fmt.Fprintf(fileTwo, "@%s\n%s\n+\n%s\n", fq.Rev.Name, dna.BasesToString(fq.Rev.Seq), QualString(fq.Rev.Qual))
	common.ExitIfError(err)
}

// WritePair takes two filenames and writes PairedEnd reads to the respective outputs.
func WritePair(readOne string, readTwo string, records []PairedEnd) {
	fileOne := fileio.EasyCreate(readOne)
	fileTwo := fileio.EasyCreate(readTwo)
	defer fileOne.Close()
	defer fileTwo.Close()
	for _, fq := range records {
		WritingHelper(fileOne, fileTwo, fq)
	}
}

// WritingChan takes two filenames and writes PairedEnd structs from an input channel and writes them to the respective files.
// Note that this function doesn't close the output channel and that this must be done when the function is called.
func WritingChan(readOne string, readTwo string, output <-chan PairedEnd, wg *sync.WaitGroup) {
	fileOne, fileTwo := fileio.EasyCreate(readOne), fileio.EasyCreate(readTwo)
	for fq := range output {
		WritingHelper(fileOne, fileTwo, fq)
	}
	fileOne.Close()
	fileTwo.Close()
	wg.Done()
}

func GoWriteFqPair(readOne string, readTwo string, data <-chan PairedEnd) {
	var wg sync.WaitGroup
	wg.Add(1)
	go WritingChan(readOne, readTwo, data, &wg)
	wg.Wait()
}

// ReadFqBigPair will take 2 readers which will convert paired end read fastq files into a paired end FastqBig struct.
// Note: while this function is return as a pointer, it's purpose is to be deferenced at the next function call.
// In addition the pointers to read one and read two will also be dereference. When sending a pointer to a struct through
// a channel, or a struct with pointers inside, memory allocated will be placed on the heap hindering performance.
func ReadFqBigPair(readerOne *fileio.ByteReader, readerTwo *fileio.ByteReader) (PairedEndBig, bool) {
	var doneOne, doneTwo bool
	var fqOne, fqTwo FastqBig = FastqBig{}, FastqBig{}
	fqOne, doneOne = ReadFqBig(readerOne)
	fqTwo, doneTwo = ReadFqBig(readerTwo)
	if !doneOne && !doneTwo {
		return PairedEndBig{Fwd: fqOne, Rev: fqTwo}, false
	} else if doneOne && doneTwo {
		return PairedEndBig{}, true
	} else {
		log.Fatalf("Error: fastq files do not end at the same time...\n")
		return PairedEndBig{}, true
	}
}
