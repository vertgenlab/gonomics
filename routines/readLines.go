package routines

import (
	"fmt"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/sam"
	"github.com/vertgenlab/gonomics/vcf"
	"io"
	"log"
	"sync"
)

type Operation func(interface{}, interface{}) interface{}

type NextLine func(reader *fileio.EasyReader) (interface{}, bool)

func goReadFile(filename string, sendLine chan interface{}, fnNext NextLine) {
	file := fileio.EasyOpen(filename)
	defer file.Close()
	var done = false
	var line interface{}

	progressMeter := 0
	for line, done = fnNext(file); done != true; line, done = fnNext(file) {
		progressMeter++
		if progressMeter % 10000 == 0 {
			log.Println("Processed", progressMeter, "Lines")
		}
		sendLine <- line
	}
	close(sendLine)
}

func GoWorkOnLine(filename string, dataStructure interface{}, fnWork Operation, fnNext NextLine, threads int) (chan interface{}, sync.WaitGroup) {
	var wg sync.WaitGroup
	records := make(chan interface{})
	output := make(chan interface{})

	go goReadFile(filename, records, fnNext)

	for k := 0; k < threads; k++ {
		wg.Add(1)
		go func() {
			for data := range records {
				a := fnWork(data, dataStructure)
				output <- a
			}
			wg.Done()
			return
		}()
	}

	go func() {
		wg.Wait()
		close(output)
	}()

	return output, wg
}

func Wait(channel chan interface{}) {
	for range channel {}
}

// Wrappers for next line functions
func NextSamLine(reader *fileio.EasyReader) (interface{}, bool) {
	_, err := reader.Peek(1)
	if err != io.EOF {sam.ReadHeader(reader)}
	return sam.NextAlignment(reader)
}

func NextVcfLine(reader *fileio.EasyReader) (interface{}, bool) {
	return vcf.NextVcf(reader)
}

func NextString(reader *fileio.EasyReader) (interface{}, bool) {
	return fileio.EasyNextLine(reader)
}

// Example Function
func goCountForReads(read interface{}, empty interface{}) interface{} {
	// Make type assertion on interface input
	data := read.(*sam.SamAln)

	// Do some work
	if sam.IsForwardRead(data) {
		return true
	} else {
		return false
	}
}
func GoCountForReads(samFilename string) {
	var count int = 0

	channel, _ := GoWorkOnLine(samFilename, nil, goCountForReads, NextSamLine, 10)

	// Listen on the output channel
	for i := range channel {
		// Declare i (output of goCountAlignedReads) as a type bool
		if !i.(bool) {
			// If read is forward, increment count
			count++
		}
	}

	// Once the count is finished, the result can be used
	// Range automatically closes when no more data is being sent, so there is no need to use the waitgroup output
	fmt.Println("There are", count, "forward reads")
}