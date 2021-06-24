//package bed provides functions for reading, writing, and manipulating Browser Extinsible Data (BED) format files.
//More information on the BED file format can be found at https://genome.ucsc.edu/FAQ/FAQformat.html#format1
package bed

import (
	"fmt"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"io"
	"log"
	"strconv"
	"strings"
	"sync"
)

// Bed stores information about genomic regions, including their location, name, score, strand, and other annotations.
type Bed struct {
	Chrom             string
	ChromStart        int
	ChromEnd          int
	Name              string
	Score             int
	Strand            Strand
	FieldsInitialized int      //number of fields that are initialized, used for smart writing.
	Annotation        []string //long form for extra fields
}

// Strand stores strand state, which can be positive, negative, or none.
type Strand byte

const (
	Positive Strand = '+'
	Negative Strand = '-'
	None     Strand = '.'
)

// String converts a bed struct to a string so it will be automatically formatted when printing with the fmt package.
func (b Bed) String() string {
	return ToString(b, b.FieldsInitialized)
}

// ToString converts a Bed struct into a BED file format string. Useful for writing to files or printing.
func ToString(bunk Bed, fields int) string {
	switch fields {
	case 3:
		return fmt.Sprintf("%s\t%d\t%d", bunk.Chrom, bunk.ChromStart, bunk.ChromEnd)
	case 4:
		return fmt.Sprintf("%s\t%d\t%d\t%s", bunk.Chrom, bunk.ChromStart, bunk.ChromEnd, bunk.Name)
	case 5:
		return fmt.Sprintf("%s\t%d\t%d\t%s\t%d", bunk.Chrom, bunk.ChromStart, bunk.ChromEnd, bunk.Name, bunk.Score)
	case 6:
		return fmt.Sprintf("%s\t%d\t%d\t%s\t%d\t%c", bunk.Chrom, bunk.ChromStart, bunk.ChromEnd, bunk.Name, bunk.Score, bunk.Strand)
	case 7:
		var out string = fmt.Sprintf("%s\t%d\t%d\t%s\t%d\t%c", bunk.Chrom, bunk.ChromStart, bunk.ChromEnd, bunk.Name, bunk.Score, bunk.Strand)
		for i := 0; i < len(bunk.Annotation); i++ {
			out = fmt.Sprintf("%s\t%s", out, bunk.Annotation[i])
		}
		return out
	default:
		log.Fatalf("Error: expecting a request to print 3 to 7 bed fields, but got: %d\n", fields)
	}
	return ""
}

// WriteBed writes an input Bed struct to an io.Writer.
func WriteBed(file io.Writer, input Bed) {
	var err error
	_, err = fmt.Fprintf(file, "%s\n", input)
	exception.PanicOnErr(err)
}

// WriteToFileHandle writes an input Bed struct to an io.Writer
func WriteToFileHandle(file io.Writer, rec Bed) {
	var err error
	_, err = fmt.Fprintf(file, "%s\n", rec)
	exception.PanicOnErr(err)
}

//Write writes a slice of Bed structs to a specified filename.
func Write(filename string, records []Bed) {
	var err error
	file := fileio.EasyCreate(filename)

	for i := range records {
		WriteToFileHandle(file, records[i])
	}
	err = file.Close()
	exception.PanicOnErr(err)
}

//Read returns a slice of Bed structs from an input filename.
func Read(filename string) []Bed {
	var line string
	var answer []Bed
	var err error
	var doneReading bool

	file := fileio.EasyOpen(filename)

	for line, doneReading = fileio.EasyNextRealLine(file); !doneReading; line, doneReading = fileio.EasyNextRealLine(file) {
		current := processBedLine(line)
		answer = append(answer, current)
	}
	err = file.Close()
	exception.PanicOnErr(err)
	return answer
}

//processBedLine is a helper function of Read that returns a Bed struct from an input line of a file.
func processBedLine(line string) Bed {
	words := strings.Split(line, "\t")
	startNum, err := strconv.Atoi(words[1])
	exception.PanicOnErr(err)
	endNum, err := strconv.Atoi(words[2])
	exception.PanicOnErr(err)

	current := Bed{Chrom: words[0], ChromStart: startNum, ChromEnd: endNum, Strand: None, FieldsInitialized: len(words)}
	if len(words) >= 4 {
		current.Name = words[3]
	}
	if len(words) >= 5 {
		current.Score, err = strconv.Atoi(words[4])
		exception.PanicOnErr(err)
	}
	if len(words) >= 6 {
		current.Strand = StringToStrand(words[5])
	}
	if len(words) >= 7 {
		for i := 6; i < len(words); i++ {
			current.Annotation = append(current.Annotation, words[i])
		}
	}
	return current
}

//StringToStrand parses a bed.Strand struct from an input string.
func StringToStrand(s string) Strand {
	switch s {
	case "+":
		return Positive
	case "-":
		return Negative
	case ".":
		return None
	default:
		log.Fatalf("Error: expected %s to be a strand that is either '+', '-', or '.'.\n", s)
		return None
	}
}

//NextBed returns a Bed struct from an input fileio.EasyReader. Returns a bool that is true when the reader is done.
func NextBed(reader *fileio.EasyReader) (Bed, bool) {
	line, done := fileio.EasyNextLine(reader)
	if done {
		return Bed{}, true
	}
	return processBedLine(line), false
}

//ReadToChan reads from a fileio.EasyReader to send Bed structs to a chan<- Bed.
func ReadToChan(file *fileio.EasyReader, data chan<- Bed, wg *sync.WaitGroup) {
	for curr, done := NextBed(file); !done; curr, done = NextBed(file) {
		data <- curr
	}
	err := file.Close()
	exception.PanicOnErr(err)
	wg.Done()
}

//GoReadToChan reads Bed entries from an input filename to a <-chan *Bed.
func GoReadToChan(filename string) <-chan Bed {
	file := fileio.EasyOpen(filename)
	var wg sync.WaitGroup
	data := make(chan Bed)
	wg.Add(1)
	go ReadToChan(file, data, &wg)

	go func() {
		wg.Wait()
		close(data)
	}()

	return data
}
