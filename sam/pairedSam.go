package sam

import (
	"fmt"
	"github.com/vertgenlab/gonomics/common"
	"os"
	"sync"
)

type PairedSamAln struct {
	FwdSam *SamAln
	RevSam *SamAln
}

func SamChanPairToFile(incomingSams <-chan *PairedSamAln, filename string, header *SamHeader, wg *sync.WaitGroup) {
	file, _ := os.Create(filename)
	defer file.Close()
	WriteHeaderToFileHandle(file, header)
	for alignedRead := range incomingSams {
		WriteAlnPairToFileHandle(file, alignedRead)
	}
	wg.Done()
}

func WriteAlnPairToFileHandle(file *os.File, aln *PairedSamAln) {
	_, err := fmt.Fprintf(file, "%s\n%s\n", SamAlnToString(aln.FwdSam), SamAlnToString(aln.RevSam))
	common.ExitIfError(err)
}
