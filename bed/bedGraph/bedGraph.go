package bedGraph

import (
	"fmt"
	"github.com/vertgenlab/gonomics/common"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"io"
	"strings"
	"sync"
)

// BedGraph stores information about genomic regions, including their location and an associated float DataValue.
// As with Bed, coordinates are 0-based, half-open.
type BedGraph struct {
	Chrom      string
	ChromStart int
	ChromEnd   int
	DataValue  float64
}

// ToString converts a BedGraph struct into a bedGraph file format string. Useful for writing to files or printing.
func ToString(b BedGraph) string {
	return fmt.Sprintf("%s\t%d\t%d\t%g", b.Chrom, b.ChromStart, b.ChromEnd, b.DataValue)
}

//Write writes a slice of BedGraph structs to a specified filename.
func Write(filename string, records []BedGraph) {
	var err error
	file := fileio.EasyCreate(filename)

	for i := range records {
		WriteToFileHandle(file, records[i])
	}
	err = file.Close()
	exception.PanicOnErr(err)
}

// WriteToFileHandle writes an input BedGraph struct to an io.Writer
func WriteToFileHandle(file io.Writer, rec BedGraph) {
	var err error
	_, err = fmt.Fprintf(file, "%s\n", rec)
	exception.PanicOnErr(err)
}

// Read returns a slice of BedGraph structs from an input filename.
func Read(filename string) []BedGraph {
	var line string
	var answer []BedGraph
	var err error
	var doneReading bool

	file := fileio.EasyOpen(filename)

	for line, doneReading = fileio.EasyNextRealLine(file); !doneReading; line, doneReading = fileio.EasyNextRealLine(file) {
		current := processBedGraphLine(line)
		answer = append(answer, current)
	}
	err = file.Close()
	exception.PanicOnErr(err)
	return answer
}

// processBedGraphLine is a helper function of Read that returns a BedGraph struct from an input string.
func processBedGraphLine(line string) BedGraph {
	words := strings.Split(line, "\t")
	startNum := common.StringToInt(words[1])
	endNum := common.StringToInt(words[2])
	dataValue := common.StringToFloat64(words[3])
	return BedGraph{Chrom: words[0], ChromStart: startNum, ChromEnd: endNum, DataValue: dataValue}
}

//NextBedGraph returns a BedGraph struct from an input fileio.EasyReader. Returns a bool that is true when the reader is done.
func NextBedGraph(reader *fileio.EasyReader) (BedGraph, bool) {
	line, done := fileio.EasyNextLine(reader)
	if done {
		return BedGraph{}, true
	}
	return processBedGraphLine(line), false
}

//ReadToChan reads from a fileio.EasyReader to send BedGraph structs to a chan<- BedGraph.
func ReadToChan(file *fileio.EasyReader, data chan<- BedGraph, wg *sync.WaitGroup) {
	for curr, done := NextBedGraph(file); !done; curr, done = NextBedGraph(file) {
		data <- curr
	}
	err := file.Close()
	exception.PanicOnErr(err)
	wg.Done()
}

//GoReadToChan reads BedGraph entries from an input filename to a <-chan BedGraph.
func GoReadToChan(filename string) <-chan BedGraph {
	file := fileio.EasyOpen(filename)
	var wg sync.WaitGroup
	data := make(chan BedGraph)
	wg.Add(1)
	go ReadToChan(file, data, &wg)

	go func() {
		wg.Wait()
		close(data)
	}()

	return data
}
