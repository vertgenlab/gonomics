//package wig provides functions to read, write, and manipulate wig files.
//more information on the WIG file format can be found at https://genome.ucsc.edu/goldenPath/help/wiggle.html

package wig

import (
	"errors"
	"fmt"
	"github.com/vertgenlab/gonomics/common"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"io"
	"log"
	"strings"
)

// Wig stores information on the chromosome location and step properties of Wig data. Individual wig values are stored in the underlying WigValue struct. Can only handle fixedStep wigs.
type Wig struct {
	StepType string
	Chrom    string
	Start    int
	Step     int
	Span     int
	Values   []float64
}

//TODO: ReadToChan() and GoReadToChan () ; write these. Check out fastq or bed or sam for examples of these ReadToChan and GoReadToChan(), these are the memory
//efficient ways to read things in. Funny things with weight groups (which handle the order of things) (concurrency in Go youtube).

// Read generates a Wig data structure from an input filename, provided as a string for a WIG format file.
func Read(filename string) []Wig {
	var curr Wig
	var finalWig []Wig
	var doneReading bool = false

	file := fileio.EasyOpen(filename)

	for curr, doneReading = NextWig(file); !doneReading; curr, doneReading = NextWig(file) { //TODO: use channels here instead of appending. Line 44 would be changed
		finalWig = append(finalWig, curr)
	}
	var err error
	err = file.Close()
	exception.PanicOnErr(err)
	return finalWig
}

func NextWig(file *fileio.EasyReader) (Wig, bool) {
	var line string
	var currentWig Wig
	var doneReading bool
	var peek []byte
	var err error
	var lineFields, chromList, startList, stepList, spanList []string

	for line, doneReading = fileio.EasyNextRealLine(file); !doneReading; line, doneReading = fileio.EasyNextRealLine(file) {
		if strings.HasPrefix(line, "fixedStep") {
			lineFields = strings.Fields(line)
			if len(lineFields) > 5 || len(lineFields) < 4 {
				log.Fatalf("Invalid number of arguments, expecting 4 or 5, received %d\n", len(lineFields))
			}
			currentWig.StepType = "fixedStep"
			chromList = strings.Split(lineFields[1], "=")
			currentWig.Chrom = chromList[1]
			startList = strings.Split(lineFields[2], "=")
			currentWig.Start = common.StringToInt(startList[1])
			stepList = strings.Split(lineFields[3], "=")
			currentWig.Step = common.StringToInt(stepList[1])
			if len(lineFields) == 5 {
				spanList = strings.Split(lineFields[4], "=")
				currentWig.Span = common.StringToInt(spanList[1])
			} else {
				currentWig.Span = -1 //signify missing
			}
		} else if strings.HasPrefix(line, "variableStep") {
			log.Fatalf("ERROR: %s is variableStep Wig, must convert to fixedStep before reading in Wig to gonomics", file.File.Name())
		} else {
			if currentWig.StepType == "" {
				log.Fatalf("ERROR: %s is missing a wig header (e.g. fixedStep chrom=chr...)", file.File.Name())
			}
			currentWig.Values = append(currentWig.Values, common.StringToFloat64(line))
		}

		peek, err = fileio.EasyPeekReal(file, 1)
		if errors.Is(err, io.EOF) {
			return currentWig, doneReading // we could return true here, but we do not in order to match the behavior of other Next functions.
		} else {
			exception.PanicOnErr(err)
		}
		if peek[0] == 'f' {
			return currentWig, doneReading
		}
	}
	return currentWig, doneReading
}

// Prints the first record in a Wig struct. Mainly used for debugging.
func PrintFirst(rec []Wig) {
	if len(rec) == 0 {
		fmt.Println("Empty Wig; length of input was zero.")
	} else {
		fmt.Printf("StepType=%s Chrom=%s Start=%d Step=%d\n", rec[0].StepType, rec[0].Chrom,
			rec[0].Start, rec[0].Step)
		if rec[0].StepType == "fixedStep" {
			fmt.Println(rec[0].Values[0])
		} else if rec[0].StepType == "variableStep" {
			log.Fatalf("ERROR: Wiggle at chrom %s and starting %d is variableStep Wig, must convert to fixedStep before reading in Wig to gonomics", rec[0].Chrom, rec[0].Start)
		}
	}
}

// Write writes a Wig data structure to a WIG format file at the input filename
func Write(filename string, rec []Wig) {
	file := fileio.EasyCreate(filename)
	for i := range rec {
		//log.Printf("Printing wig object: %d\n", i)
		WriteToFileHandle(file, rec[i])
	}
	var err error
	err = file.Close()
	exception.PanicOnErr(err)
}

// WriteToFileHandle is an helper function for Write that writes the Wig data structure to a file.
func WriteToFileHandle(file io.Writer, rec Wig) {
	var err error
	if rec.StepType == "fixedStep" {
		if rec.Span != -1 { //if there was a span in the header line
			_, err = fmt.Fprintf(file, "%s chrom=%s start=%d step=%d span=%d\n", rec.StepType, rec.Chrom,
				rec.Start, rec.Step, rec.Span)
		} else {
			_, err = fmt.Fprintf(file, "%s chrom=%s start=%d step=%d\n", rec.StepType, rec.Chrom,
				rec.Start, rec.Step)
		}
		common.ExitIfError(err)
	} else if rec.StepType == "variableStep" {
		log.Fatalf("ERROR: %s is variableStep Wig, must convert to fixedStep before reading in Wig to gonomics", rec.StepType)
	} else {
		log.Fatalf("Invalid step type for wig.")
	}

	for i := range rec.Values {
		if rec.StepType == "fixedStep" {
			if rec.Values[i] == 0 {
				_, err = fmt.Fprintf(file, "0\n") //will turn a 0.000000 to a 0 to save mem
			} else {
				_, err = fmt.Fprintf(file, "%f\n", rec.Values[i])
			} //Only print significant figures to some capacity? rpkm
			// We want to make this as concise as possible.
			// How can we break up a wig into sections that actually has data and skip over large sections of zeros?
			common.ExitIfError(err)
		}

	}
}

// ChromToSlice returns the values from a wig entry matching a user-specified chromosome name
func ChromToSlice(w []Wig, chrom string) []float64 {
	var output []float64
	for _, v := range w {
		if v.Chrom == chrom {
			if v.Step == 1 {
				if v.Start == 1 {
					output = make([]float64, len(v.Values))
					output = v.Values
				} else {
					log.Fatalf("Invalid start coordinate.")
				}
			} else {
				log.Fatalf("invalid step size, step size must be 1.")
			}
		}
	}
	return output
}
