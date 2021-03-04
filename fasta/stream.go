package fasta

import (
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"sync"
)

// ReadToChan is a helper function of GoReadToChan.
func ReadToChan(file *fileio.EasyReader, data chan<- Fasta, wg *sync.WaitGroup) {
	for curr, done := NextFasta(file); !done; curr, done = NextFasta(file) {
		data <- curr
	}
	exception.PanicOnErr(file.Close())
	wg.Done()
}

// GoReadToChan reads fasta records from an input filename and returns a channel of Fasta structs.
func GoReadToChan(filename string) <-chan Fasta {
	file := fileio.EasyOpen(filename)
	var wg sync.WaitGroup
	data := make(chan Fasta)
	wg.Add(1)
	go ReadToChan(file, data, &wg)

	go func() {
		wg.Wait()
		close(data)
	}()

	return data
}
