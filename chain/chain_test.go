package chain

import (
	"github.com/vertgenlab/gonomics/fileio"
	"log"
	"os"
	"strings"
	"sync"
	"testing"
)

var dir []string = []string{"testdata/big.chain", "testdata/twoChainZ.chain"}

func TestReaderAndWriter(t *testing.T) {
	for _, readFile := range dir {
		//1st pass, test standard read and write (non-channel) functuions
		writeToFile := readFile + ".tmp"
		readerChains, hashTags := Read(readFile)
		Write(writeToFile, readerChains, hashTags)

		//2nd pass: read back and write using goroutines
		file := fileio.EasyOpen(writeToFile)
		reader := make(chan *Chain)
		defer file.Close()
		goHashTags := ReadHeaderComments(file)
		go ReadToChan(file, reader)
		testGoRoutines := "testdata/goReaderWriter.chain"
		var wg sync.WaitGroup
		wg.Add(1)
		go WriteToFile(testGoRoutines, reader, goHashTags, &wg)
		wg.Wait()

		//Read back in file that was processeds with goroutines using the standard version. Could read back in using goroutines, but this was faster to set-up
		testData, _ := Read(testGoRoutines)
		if !Equal(readerChains, testData) {
			t.Errorf("Error: File read in did not match the file writing out")
		} else {
			log.Printf("%s: Passed simple read and write test!\n", strings.Split(readFile, "/")[1])
		}
		os.Remove(writeToFile)
		os.Remove(testGoRoutines)
	}
}
