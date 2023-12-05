// Package wig provides functions to read, write, and manipulate wig files.
// more information on the WIG file format can be found at https://genome.ucsc.edu/goldenPath/help/wiggle.html
package wig

import (
	"errors"
	"fmt"
	"github.com/vertgenlab/gonomics/chromInfo"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/numbers/parse"
	"io"
	"log"
	"sort"
	"strings"
	"sync"
)

// Wig stores information on the chromosome location and step properties of Wig data. Individual wig values
// are stored in the underlying WigValue struct. Can only handle fixedStep wigs.
type Wig struct {
	StepType string
	Chrom    string
	Start    int
	Step     int
	Span     int
	Values   []float64
}

// NextWig returns a Wig struct from an input fileio.EasyReader. Returns a bool that is true when the reader is done.
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
			currentWig.Start = parse.StringToInt(startList[1])
			stepList = strings.Split(lineFields[3], "=")
			currentWig.Step = parse.StringToInt(stepList[1])
			if len(lineFields) == 5 {
				spanList = strings.Split(lineFields[4], "=")
				currentWig.Span = parse.StringToInt(spanList[1])
			} else {
				currentWig.Span = -1 //signify missing
			}
		} else if strings.HasPrefix(line, "variableStep") {
			log.Fatalf("ERROR: %s is variableStep Wig, must convert to fixedStep before reading in Wig to gonomics", file.File.Name())
		} else {
			if currentWig.StepType == "" {
				log.Fatalf("ERROR: %s is missing a wig header (e.g. fixedStep chrom=chr...)", file.File.Name())
			}
			currentWig.Values = append(currentWig.Values, parse.StringToFloat64(line))
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

// ReadToChan reads from a fileio.EasyReader to send Wig structs to a chan<- Wig.
// When it has finished reading the file, ReadToChan will call Done on the waitgroup.
func ReadToChan(file *fileio.EasyReader, data chan<- Wig, wg *sync.WaitGroup) {
	for curr, done := NextWig(file); !done; curr, done = NextWig(file) {
		data <- curr
	}
	err := file.Close()
	exception.PanicOnErr(err)
	wg.Done()
}

// GoReadToChan reads Wig entries from an input filename to a <-chan Wig.
func GoReadToChan(filename string) <-chan Wig {
	file := fileio.EasyOpen(filename)
	var wg sync.WaitGroup
	data := make(chan Wig)
	wg.Add(1)
	go ReadToChan(file, data, &wg)

	go func() {
		wg.Wait()
		close(data)
	}()

	return data
}

// WriteMap writes an input map[string]Wig to an output filename. Entries in the map
// will be written in alphanumeric order of the map keys.
func WriteMap(filename string, records map[string]Wig) {
	var err error
	keys := make([]string, 0)
	out := fileio.EasyCreate(filename)
	for currKey, _ := range records {
		keys = append(keys, currKey)
	}
	//sorting by keys to write in deterministic order
	sort.Strings(keys)
	for _, currKey := range keys {
		WriteToFileHandle(out, records[currKey])
	}
	err = out.Close()
	exception.PanicOnErr(err)
}

// WriteToFileHandle is an helper function for Write that writes the Wig data structure to an io.Writer
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
		exception.PanicOnErr(err)
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
			exception.PanicOnErr(err)
		}
	}
}

// ReadWholeGenome creates a whole genome wig map, where each key is a chromosome name and the value is the corresponding
// wig struct. The wig struct values slice is the size of the whole chromosome, as specified by an input chromSizeFile.
// Positions in the genome where there is no wig coverage are set to a user-specified defaultValue.
// This function is robust against different step sizes between entries and multiple wig entries per chromosome, but
// is currently limited to 'fixedStep' beds.
func ReadWholeGenome(filename string, chromSizeFile string, defaultValue float64) map[string]Wig {
	var currValue, currPos, currOffset int
	var foundInMap bool
	sizes := chromInfo.ReadToMap(chromSizeFile)
	answer := MakeSkeleton(sizes, defaultValue)
	wigChan := GoReadToChan(filename)
	for currWig := range wigChan {
		if _, foundInMap = answer[currWig.Chrom]; !foundInMap {
			log.Fatalf("Error: chrom name in wig file: %s, not found in reference genome chrom sizes.\n", currWig.Chrom)
		}
		currPos = currWig.Start - 1 // from 1-based to zero-based
		for currValue = range currWig.Values {
			if currPos > len(answer[currWig.Chrom].Values) {
				log.Fatalf("Error: position of values in the input wig exceed the chrom length specified in the chrom sizes file. Offending entry on Chr: %v at start: %v\n", currWig.Chrom, currWig.Start)
			}
			for currOffset = 0; currOffset < currWig.Step; currOffset++ {
				//TODO: Check if value is not defaultValue first and fatal
				answer[currWig.Chrom].Values[currPos] = currWig.Values[currValue]
				currPos++
			}
		}
	}
	return answer
}

// MakeSkeleton creates a fixed-step map[string]Wig, with one Wig entry per chromosome, where the size of the
// Values slice in each struct is equal to the chromosome length, as specified by an input chromSizes map[string]chromInfo.ChromInfo.
// 'Skeleton' refers to the fact that the values are uniformly set to a defaultValue, and contain no information.
func MakeSkeleton(chromSizes map[string]chromInfo.ChromInfo, defaultValue float64) map[string]Wig {
	var answer = make(map[string]Wig, 0)
	var currentWig Wig
	var x int

	for _, v := range chromSizes {
		currentWig = Wig{StepType: "fixedStep", Chrom: v.Name, Start: 1, Step: 1, Span: -1}
		currentWig.Values = make([]float64, v.Size)
		for x = 0; x < v.Size; x++ {
			currentWig.Values[x] = defaultValue
		}
		answer[v.Name] = currentWig
	}
	return answer
}
