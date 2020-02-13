package sam

import (
	"fmt"
	"github.com/vertgenlab/gonomics/common"
	"log"
	"os"
	"sync"
)

type PairedSamAln struct {
	FwdSam *SamAln
	RevSam *SamAln
}

func SamChanPairToFile(incomingSams <-chan *PairedSamAln, file *os.File, wg *sync.WaitGroup) {

	for alignedRead := range incomingSams {
		WriteAlnPairToFileHandle(file, alignedRead)
	}
	log.Printf("Finished aligning read pair!! %s\n")
	wg.Done()
}

func WriteAlnPairToFileHandle(file *os.File, aln *PairedSamAln) {
	_, err := fmt.Fprintf(file, "%s\n%s\n", SamAlnToString(aln.FwdSam), SamAlnToString(aln.RevSam))

	common.ExitIfError(err)
}
