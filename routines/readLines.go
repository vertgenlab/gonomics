package routines

import (
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

func GoWorkOnLine(filename string, dataStructure interface{}, fnWork Operation, fnNext NextLine, threads int) chan interface{} {
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

	return output
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