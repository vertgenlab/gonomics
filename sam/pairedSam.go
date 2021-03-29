package sam

import (
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"sync"
)

// PairedSamAln wraps to mate pair alignments into a single struct.
type PairedSamAln struct {
	FwdSam Aln
	RevSam Aln
}

// SamChaPairToFile writes a channel of PairedSamAln to a file.
func SamChanPairToFile(incomingSams <-chan PairedSamAln, filename string, header Header, wg *sync.WaitGroup) {
	file := fileio.EasyCreate(filename)
	if header.Text != nil {
		WriteHeaderToFileHandle(file, header)
	}
	for alignedRead := range incomingSams {
		WriteToFileHandle(file, alignedRead.FwdSam)
		WriteToFileHandle(file, alignedRead.RevSam)
	}
	err := file.Close()
	exception.PanicOnErr(err)
	wg.Done()
}