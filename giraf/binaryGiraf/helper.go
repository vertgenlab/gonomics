package binaryGiraf

import (
	"io"
	"sync"

	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/genomeGraph"
	"github.com/vertgenlab/gonomics/giraf"
)

// Read will input a giraf.fe file and will decompress into a slice of giraf.Giraf.
// Requires the graph (.gg) file used to generate the giraf for decompression.
func Read(giraffeFile string, graphFile string) []giraf.Giraf {
	var answer []giraf.Giraf
	graph := genomeGraph.Read(graphFile)
	reader := NewBinReader(giraffeFile)
	var curr giraf.Giraf
	var err error
	for curr, err = ReadGiraf(reader, graph); err == nil; curr, err = ReadGiraf(reader, graph) {
		answer = append(answer, curr)
	}

	if err != io.EOF {
		exception.PanicOnErr(err)
	}

	// Close reader
	err = reader.bg.Close()
	exception.PanicOnErr(err)

	return answer
}

// ReadToChan will input a fileio.EasyReader and decompress girafs in the file into a stream (channel) of giraf.Giraf.
// ReadToChan is designed to be run as a goroutine. See GoReadToChan for internally handled goroutines.
func ReadToChan(file *fileio.EasyReader, graph *genomeGraph.GenomeGraph, data chan<- giraf.Giraf, wg *sync.WaitGroup) {
	var curr giraf.Giraf
	var err error
	reader := NewBinReader(file.File.Name())
	for curr, err = ReadGiraf(reader, graph); err == nil; curr, err = ReadGiraf(reader, graph) {
		data <- curr
	}

	if err != io.EOF {
		exception.PanicOnErr(err)
	}

	// Close reader
	err = reader.bg.Close()
	exception.PanicOnErr(err)

	err = file.Close()
	exception.PanicOnErr(err)
	wg.Done()
}

// GoReadToChan is a wrapper for the ReadToChan function that handles goroutine creation internally. Decompresses input giraf.fe file
// into a buffered channel (len 1024) of giraf.Giraf records.
func GoReadToChan(giraffeFile string, graphFile string) <-chan giraf.Giraf {
	file := fileio.EasyOpen(giraffeFile)
	graph := genomeGraph.Read(graphFile)
	var wg sync.WaitGroup
	data := make(chan giraf.Giraf, 1024)
	wg.Add(1)
	go ReadToChan(file, graph, data, &wg)

	go func() {
		wg.Wait()
		close(data)
	}()

	return data
}

// GirafChanToBinary inputs a channel of *giraf.Giraf and compresses records from the stream into a output giraf.fe file.
// TODO: move GoReadToChan functions in the giraf package from *Giraf to Giraf.
func GirafChanToBinary(filename string, input <-chan *giraf.Giraf, wg *sync.WaitGroup) {
	file := fileio.EasyCreate(filename)
	writer := NewBinWriter(file)
	var err error

	for record := range input {
		err = WriteGiraf(writer, record)
		exception.PanicOnErr(err)
	}
	err = writer.bg.Close()
	exception.PanicOnErr(err)
	wg.Done()

	err = file.Close()
	exception.PanicOnErr(err)
}
