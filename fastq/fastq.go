// Package fastq provides functions for reading, writing, and manipulating data in the fastq file format.
package fastq

import (
	"fmt"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"log"
	"sync"
)

// Fastq stores the name, sequence, and quality scores for each 4 line record of a Fastq file.
type Fastq struct {
	Name string
	Seq  []dna.Base
	Qual []uint8
}

// Read sends all records in a fastq format file to a []Fastq.
func Read(filename string) []Fastq {
	file := fileio.EasyOpen(filename)
	answer := ReadFastqs(file)
	err := file.Close()
	exception.PanicOnErr(err)
	return answer
}

// ReadToChan reads from an EasyReader, sending fastq records to data.  When it gest to the end
// of the fastq file, it will call Done() on the waitgroup.
func ReadToChan(file *fileio.EasyReader, data chan<- Fastq, wg *sync.WaitGroup) {
	for curr, done := NextFastq(file); !done; curr, done = NextFastq(file) {
		data <- curr
	}
	err := file.Close()
	exception.PanicOnErr(err)
	wg.Done()
}

// GoReadToChan launches a go routine that parses Fastq structs from an input filename and returns a chan of those Fastq structs.
func GoReadToChan(filename string) <-chan Fastq {
	file := fileio.EasyOpen(filename)
	var wg sync.WaitGroup
	data := make(chan Fastq)
	wg.Add(1)
	go ReadToChan(file, data, &wg)

	go func() {
		wg.Wait()
		close(data)
	}()

	return data
}

// Write writes all records in an input []Fastq to a file at an input filename string.
func Write(filename string, records []Fastq) {
	file := fileio.EasyCreate(filename)
	for _, fq := range records {
		WriteToFileHandle(file, fq)
	}
	err := file.Close()
	exception.PanicOnErr(err)
}

// WriteToFileHandle writes an individual fastq record to a given fileio.EasyWriter.
func WriteToFileHandle(file *fileio.EasyWriter, fq Fastq) {
	var err error
	_, err = fmt.Fprintf(file, "@%s\n%s\n+\n%s\n", fq.Name, dna.BasesToString(fq.Seq), QualString(fq.Qual))
	exception.PanicOnErr(err)
}

// processFastqRecord is a helper function of NextFastq that parses a fastq struct from an input line of a fastq file provided as a string.
func processFastqRecord(line1 string, line2 string, line3 string, line4 string) Fastq {
	var curr Fastq
	if line3 != "+" {
		log.Fatalf("Error: This line should be a + (plus) sign \n")
	}
	curr = Fastq{Name: line1[1:], Seq: dna.StringToBases(line2), Qual: ToQual([]byte(line4))}
	return curr
}

// NextFastq is a helper function of Read and GoReadToChan. NextFastq checks a reader for additional lines of data and returns a fastq struct
// and a bool telling the caller if it is done reading from the file.
// When it is done reading from the EasyReader, it will return an empty fastq and true.
func NextFastq(reader *fileio.EasyReader) (Fastq, bool) {
	line, done := fileio.EasyNextLine(reader)
	line2, done2 := fileio.EasyNextLine(reader)
	line3, done3 := fileio.EasyNextLine(reader)
	line4, done4 := fileio.EasyNextLine(reader)
	if done {
		return Fastq{}, true
	}
	if done2 || done3 || done4 {
		log.Fatalf("Error: There is an empty line in this fastq record\n")
	}
	return processFastqRecord(line, line2, line3, line4), false
}

// ReadFastqs reads a slice of fastq from a fileio.EasyReader.
func ReadFastqs(er *fileio.EasyReader) []Fastq {
	var curr Fastq
	var done bool
	var answer []Fastq
	for curr, done = NextFastq(er); !done; curr, done = NextFastq(er) {
		answer = append(answer, curr)
	}
	return answer
}

// String converts a fastq struct to a string so it will be automatically formatted when printing with the fmt format.
func (fq *Fastq) String() string {
	return fmt.Sprintf("@%s\n%s\n+\n%s\n", fq.Name, dna.BasesToString(fq.Seq), QualString(fq.Qual))
}
