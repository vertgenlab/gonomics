// Package bedpe provides functions for handling Bed Paired End format files, used for
// handling pairs of genomic regions.

package bedpe

import (
	"fmt"
	"io"
	"log"
	"strings"
	"sync"

	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/common"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
)

type BedPe struct {
	A bed.Bed
	B bed.Bed
}

type BedPeHalf struct {
	Chrom      string
	ChromStart int
	ChromEnd   int
	Home       *BedPe
}

// String converts a BedPe struct to a string so it will be automatically formatted when printing with the fmt package.
func (b BedPe) String() string {
	return ToString(b, b.A.FieldsInitialized)
}

// ToString converts a BedPe struct into a BedPe file format string. Useful for writing to files or printing.
func ToString(bunk BedPe, fields int) string {
	switch {
	case fields == 6:
		return fmt.Sprintf("%s\t%d\t%d\t%s\t%d\t%d", bunk.A.Chrom, bunk.A.ChromStart, bunk.A.ChromEnd, bunk.B.Chrom, bunk.B.ChromStart, bunk.B.ChromEnd)
	case fields == 7:
		return fmt.Sprintf("%s\t%d\t%d\t%s\t%d\t%d\t%s", bunk.A.Chrom, bunk.A.ChromStart, bunk.A.ChromEnd, bunk.B.Chrom, bunk.B.ChromStart, bunk.B.ChromEnd, bunk.A.Name)
	case fields == 8:
		return fmt.Sprintf("%s\t%d\t%d\t%s\t%d\t%d\t%s\t%d", bunk.A.Chrom, bunk.A.ChromStart, bunk.A.ChromEnd, bunk.B.Chrom, bunk.B.ChromStart, bunk.B.ChromEnd, bunk.A.Name, bunk.A.Score)
	case fields == 9:
		return fmt.Sprintf("%s\t%d\t%d\t%s\t%d\t%d\t%s\t%d\t%c", bunk.A.Chrom, bunk.A.ChromStart, bunk.A.ChromEnd, bunk.B.Chrom, bunk.B.ChromStart, bunk.B.ChromEnd, bunk.A.Name, bunk.A.Score, bunk.A.Strand)
	case fields == 10:
		return fmt.Sprintf("%s\t%d\t%d\t%s\t%d\t%d\t%s\t%d\t%c\t%c", bunk.A.Chrom, bunk.A.ChromStart, bunk.A.ChromEnd, bunk.B.Chrom, bunk.B.ChromStart, bunk.B.ChromEnd, bunk.A.Name, bunk.A.Score, bunk.A.Strand, bunk.B.Strand)
	case fields >= 11:
		var out string = fmt.Sprintf("%s\t%d\t%d\t%s\t%d\t%d\t%s\t%d\t%c\t%c", bunk.A.Chrom, bunk.A.ChromStart, bunk.A.ChromEnd, bunk.B.Chrom, bunk.B.ChromStart, bunk.B.ChromEnd, bunk.A.Name, bunk.A.Score, bunk.A.Strand, bunk.B.Strand)
		for i := range bunk.A.Annotation {
			out = fmt.Sprintf("%s\t%s", out, bunk.A.Annotation[i])
		}
		return out
	default:
		log.Fatalf("Error: expecting a request to print at least 6 bedPe fields, but got: %d\n", fields)
	}
	return ""
}

// WriteToFileHandle writes an input BedPe struct to an io.Writer.
func WriteToFileHandle(file io.Writer, rec BedPe) {
	var err error
	_, err = fmt.Fprintf(file, "%s\n", rec)
	exception.PanicOnErr(err)
}

// Write writes a slice of BedPe structs to a specified filename.
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
	var answer []BedPe
	var err error

	file := fileio.EasyOpen(filename)

	for curr, done := NextBedPe(file); !done; curr, done = NextBedPe(file) {
		answer = append(answer, curr)
	}
	err = file.Close()
	exception.PanicOnErr(err)
	return answer
}

// processBedPeLine is a helper function of Read that returns a BedPe struct from an input line of a file.
func processBedPeLine(line string) BedPe {
	var startANum, endANum, startBNum, endBNum int
	words := strings.Split(line, "\t")
	startANum = common.StringToInt(words[1])
	endANum = common.StringToInt(words[2])
	startBNum = common.StringToInt(words[4])
	endBNum = common.StringToInt(words[5])

	current := BedPe{A: bed.Bed{
		Chrom:             words[0],
		ChromStart:        startANum,
		ChromEnd:          endANum,
		Strand:            bed.None,
		FieldsInitialized: len(words),
	},
		B: bed.Bed{
			Chrom:             words[3],
			ChromStart:        startBNum,
			ChromEnd:          endBNum,
			Strand:            bed.None,
			FieldsInitialized: len(words),
		},
	}
	if len(words) >= 7 {
		current.A.Name, current.B.Name = words[6], words[6]
	}
	if len(words) >= 8 {
		current.A.Score, current.B.Score = common.StringToInt(words[7]), common.StringToInt(words[7])
	}
	if len(words) >= 9 {
		current.A.Strand = bed.StringToStrand(words[8])
	}
	if len(words) >= 10 {
		current.B.Strand = bed.StringToStrand(words[9])
	}
	if len(words) >= 11 {
		for i := 10; i < len(words); i++ {
			current.A.Annotation = append(current.A.Annotation, words[i])
			current.B.Annotation = append(current.B.Annotation, words[i])
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

// ReadToChan reads from a fileio.EasyReader to send BedPe structs to a chan<- BedPe.
func ReadToChan(file *fileio.EasyReader, data chan<- BedPe, wg *sync.WaitGroup) {
	for curr, done := NextBedPe(file); !done; curr, done = NextBedPe(file) {
		data <- curr
	}
	err := file.Close()
	exception.PanicOnErr(err)
	wg.Done()
}

// GoReadToChan reads BedPe entries from an input filename to a <-chan BedPe.
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

// SplitBedPe takes in a bedPe and creates two half based on A and B values in bedPe.
func SplitBedPe(in BedPe) (BedPeHalf, BedPeHalf) {
	left := BedPeHalf{
		Chrom:      in.A.Chrom,
		ChromStart: in.A.ChromStart,
		ChromEnd:   in.A.ChromEnd,
		Home:       &in,
	}
	right := BedPeHalf{
		Chrom:      in.B.Chrom,
		ChromStart: in.B.ChromStart,
		ChromEnd:   in.B.ChromEnd,
		Home:       &in,
	}

	return left, right
}
