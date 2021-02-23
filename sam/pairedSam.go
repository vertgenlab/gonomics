package sam

import (
	"fmt"
	"github.com/vertgenlab/gonomics/common"
	"github.com/vertgenlab/gonomics/fileio"
	"sync"
)

type PairedSamAln struct {
	FwdSam *SamAln
	RevSam *SamAln
}

func SamChanPairToFile(incomingSams <-chan *PairedSamAln, filename string, header *SamHeader, wg *sync.WaitGroup) {
	file := fileio.EasyCreate(filename)
	defer file.Close()
	if header != nil {
		WriteHeaderToFileHandle(file, header)
	}
	for alignedRead := range incomingSams {
		WriteAlnPairToFileHandle(file, alignedRead)
	}
	wg.Done()
}

func WriteAlnPairToFileHandle(file *fileio.EasyWriter, aln *PairedSamAln) {
	_, err := fmt.Fprintf(file, "%s\n%s\n", SamAlnToString(aln.FwdSam), SamAlnToString(aln.RevSam))
	common.ExitIfError(err)
}
