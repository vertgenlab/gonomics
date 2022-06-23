// Package bedpe provides functions for handling Bed Paired End format files, used for
// handling pairs of genomic regions.

package bedpe

import (
	"fmt"
	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/common"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"io"
	"log"
	"strings"
	"sync"
)

type BedPe struct {
	ChromA            string
	ChromStartA       int
	ChromEndA         int
	ChromB            string
	ChromStartB       int
	ChromEndB         int
	Name              string
	Score             int
	StrandA           bed.Strand
	StrandB           bed.Strand
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
		return fmt.Sprintf("%s\t%d\t%d\t%s\t%d\t%d", bunk.ChromA, bunk.ChromStartA, bunk.ChromEndA, bunk.ChromB, bunk.ChromStartB, bunk.ChromEndB)
	case fields == 7:
		return fmt.Sprintf("%s\t%d\t%d\t%s\t%d\t%d\t%s", bunk.ChromA, bunk.ChromStartA, bunk.ChromEndA, bunk.ChromB, bunk.ChromStartB, bunk.ChromEndB, bunk.Name)
	case fields == 8:
		return fmt.Sprintf("%s\t%d\t%d\t%s\t%d\t%d\t%s\t%d", bunk.ChromA, bunk.ChromStartA, bunk.ChromEndA, bunk.ChromB, bunk.ChromStartB, bunk.ChromEndB, bunk.Name, bunk.Score)
	case fields == 9:
		return fmt.Sprintf("%s\t%d\t%d\t%s\t%d\t%d\t%s\t%d\t%c", bunk.ChromA, bunk.ChromStartA, bunk.ChromEndA, bunk.ChromB, bunk.ChromStartB, bunk.ChromEndB, bunk.Name, bunk.Score, bunk.StrandA)
	case fields == 10:
		return fmt.Sprintf("%s\t%d\t%d\t%s\t%d\t%d\t%s\t%d\t%c\t%c", bunk.ChromA, bunk.ChromStartA, bunk.ChromEndA, bunk.ChromB, bunk.ChromStartB, bunk.ChromEndB, bunk.Name, bunk.Score, bunk.StrandA, bunk.StrandB)
	case fields >= 11:
		var out string = fmt.Sprintf("%s\t%d\t%d\t%s\t%d\t%d\t%s\t%d\t%c\t%c", bunk.ChromA, bunk.ChromStartA, bunk.ChromEndA, bunk.ChromB, bunk.ChromStartB, bunk.ChromEndB, bunk.Name, bunk.Score, bunk.StrandA, bunk.StrandB)
		for i := range bunk.Annotation {
			out = fmt.Sprintf("%s\t%s", out, bunk.Annotation[i])
		}
		return out
	default:
		log.Fatalf("Error: expecting a request to print at least 6 bedPe fields, but got: %d\n", fields)
	}
	return ""
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
	var start1Num, end1Num, start2Num, end2Num int
	words := strings.Split(line, "\t")
	start1Num = common.StringToInt(words[1])
	end1Num = common.StringToInt(words[2])
	start2Num = common.StringToInt(words[4])
	end2Num = common.StringToInt(words[5])

	current := BedPe{ChromA: words[0],
		ChromStartA:       start1Num,
		ChromEndA:         end1Num,
		ChromB:            words[3],
		ChromStartB:       start2Num,
		ChromEndB:         end2Num,
		StrandA:           bed.None,
		StrandB:           bed.None,
		FieldsInitialized: len(words),
	}
	if len(words) >= 7 {
		current.Name = words[6]
	}
	if len(words) >= 8 {
		current.Score = common.StringToInt(words[7])
	}
	if len(words) >= 9 {
		current.StrandA = bed.StringToStrand(words[8])
	}
	if len(words) >= 10 {
		current.StrandB = bed.StringToStrand(words[9])
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
