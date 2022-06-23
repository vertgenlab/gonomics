// Package bedpe provides functions for handling Bed Paired End format files, used for
// handling pairs of genomic regions.

package bedpe

import (
	"fmt"
	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"io"
	"log"
	"strconv"
	"strings"
	"sync"
)

type BedPe struct {
	Chrom1            string
	ChromStart1       int
	ChromEnd1         int
	Chrom2            string
	ChromStart2       int
	ChromEnd2         int
	Name              string
	Score             int
	Strand1           bed.Strand
	Strand2           bed.Strand
	FieldsInitialized int      //number of fields that are initialized, used for smart writing.
	Annotation        []string //long form for extra fields
}

// String converts a BedPe struct to a string so it will be automatically formatted when printing with the fmt package.
func (b BedPe) String() string {
	return ToString(b, b.FieldsInitialized)
}

// ToString converts a BedPe struct into a BedPe file format string. Useful for writing to files or printing.
func ToString(bunk BedPe, fields int) string {
	switch {
	case fields == 6:
		return fmt.Sprintf("%s\t%d\t%d\t%s\t%d\t%d", bunk.Chrom1, bunk.ChromStart1, bunk.ChromEnd1, bunk.Chrom2, bunk.ChromStart2, bunk.ChromEnd2)
	case fields == 7:
		return fmt.Sprintf("%s\t%d\t%d\t%s\t%d\t%d\t%s", bunk.Chrom1, bunk.ChromStart1, bunk.ChromEnd1, bunk.Chrom2, bunk.ChromStart2, bunk.ChromEnd2, bunk.Name)
	case fields == 8:
		return fmt.Sprintf("%s\t%d\t%d\t%s\t%d\t%d\t%s\t%d", bunk.Chrom1, bunk.ChromStart1, bunk.ChromEnd1, bunk.Chrom2, bunk.ChromStart2, bunk.ChromEnd2, bunk.Name, bunk.Score)
	case fields == 9:
		return fmt.Sprintf("%s\t%d\t%d\t%s\t%d\t%d\t%s\t%d\t%c", bunk.Chrom1, bunk.ChromStart1, bunk.ChromEnd1, bunk.Chrom2, bunk.ChromStart2, bunk.ChromEnd2, bunk.Name, bunk.Score, bunk.Strand1)
	case fields == 10:
		return fmt.Sprintf("%s\t%d\t%d\t%s\t%d\t%d\t%s\t%d\t%c\t%c", bunk.Chrom1, bunk.ChromStart1, bunk.ChromEnd1, bunk.Chrom2, bunk.ChromStart2, bunk.ChromEnd2, bunk.Name, bunk.Score, bunk.Strand1, bunk.Strand2)
	case fields >= 11:
		var out string = fmt.Sprintf("%s\t%d\t%d\t%s\t%d\t%d\t%s\t%d\t%c\t%c", bunk.Chrom1, bunk.ChromStart1, bunk.ChromEnd1, bunk.Chrom2, bunk.ChromStart2, bunk.ChromEnd2, bunk.Name, bunk.Score, bunk.Strand1, bunk.Strand2)
		for i := range bunk.Annotation {
			out = fmt.Sprintf("%s\t%s", out, bunk.Annotation[i])
		}
		return out
	default:
		log.Fatalf("Error: expecting a request to print at least 6 bedPe fields, but got: %d\n", fields)
	}
	return ""
}

// WriteBed writes an input BedPe struct to an io.Writer.
func WriteBedPe(file io.Writer, input BedPe) {
	var err error
	_, err = fmt.Fprintf(file, "%s\n", input)
	exception.PanicOnErr(err)
}

// WriteToFileHandle writes an input BedPe struct to an io.Writer
func WriteToFileHandle(file io.Writer, rec BedPe) {
	var err error
	_, err = fmt.Fprintf(file, "%s\n", rec)
	exception.PanicOnErr(err)
}

//Write writes a slice of BedPe structs to a specified filename.
func Write(filename string, records []BedPe) {
	var err error
	file := fileio.EasyCreate(filename)

	for i := range records {
		WriteToFileHandle(file, records[i])
	}
	err = file.Close()
	exception.PanicOnErr(err)
}

// Read returns a slice of BedPe structs from an input filename.
func Read(filename string) []BedPe {
	var line string
	var answer []BedPe
	var err error
	var doneReading bool

	file := fileio.EasyOpen(filename)

	for line, doneReading = fileio.EasyNextRealLine(file); !doneReading; line, doneReading = fileio.EasyNextRealLine(file) {
		current := processBedPeLine(line)
		answer = append(answer, current)
	}
	err = file.Close()
	exception.PanicOnErr(err)
	return answer
}

//processBedPeLine is a helper function of Read that returns a BedPe struct from an input line of a file.
func processBedPeLine(line string) BedPe {
	var err error
	var start1Num, end1Num, start2Num, end2Num int
	words := strings.Split(line, "\t")
	start1Num, err = strconv.Atoi(words[1])
	exception.PanicOnErr(err)
	end1Num, err = strconv.Atoi(words[2])
	exception.PanicOnErr(err)
	start2Num, err = strconv.Atoi(words[4])
	exception.PanicOnErr(err)
	end2Num, err = strconv.Atoi(words[5])
	exception.PanicOnErr(err)
	current := BedPe{Chrom1: words[0],
		ChromStart1:       start1Num,
		ChromEnd1:         end1Num,
		Chrom2:            words[3],
		ChromStart2:       start2Num,
		ChromEnd2:         end2Num,
		Strand1:           bed.None,
		Strand2:           bed.None,
		FieldsInitialized: len(words),
	}
	if len(words) >= 7 {
		current.Name = words[6]
	}
	if len(words) >= 8 {
		current.Score, err = strconv.Atoi(words[7])
		exception.PanicOnErr(err)
	}
	if len(words) >= 9 {
		current.Strand1 = bed.StringToStrand(words[8])
	}
	if len(words) >= 10 {
		current.Strand2 = bed.StringToStrand(words[9])
	}
	if len(words) >= 11 {
		for i := 10; i < len(words); i++ {
			current.Annotation = append(current.Annotation, words[i])
		}
	}

	return current
}

// NextBedPe returns a BedPe struct from an input fileio.EasyReader. Returns a bool that is true when the reader is done.
func NextBedPe(reader *fileio.EasyReader) (BedPe, bool) {
	line, done := fileio.EasyNextLine(reader)
	if done {
		return BedPe{}, true
	}
	return processBedPeLine(line), false
}

//ReadToChan reads from a fileio.EasyReader to send BedPe structs to a chan<- BedPe.
func ReadToChan(file *fileio.EasyReader, data chan<- BedPe, wg *sync.WaitGroup) {
	for curr, done := NextBedPe(file); !done; curr, done = NextBedPe(file) {
		data <- curr
	}
	err := file.Close()
	exception.PanicOnErr(err)
	wg.Done()
}

//GoReadToChan reads BedPe entries from an input filename to a <-chan BedPe.
func GoReadToChan(filename string) <-chan BedPe {
	file := fileio.EasyOpen(filename)
	var wg sync.WaitGroup
	data := make(chan BedPe)
	wg.Add(1)
	go ReadToChan(file, data, &wg)

	go func() {
		wg.Wait()
		close(data)
	}()

	return data
}
