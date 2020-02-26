package simpleGraph

import (
	"github.com/vertgenlab/gonomics/fastq"
	"github.com/vertgenlab/gonomics/sam"
	"sync"
)

func gswWorker(gg *SimpleGraph, seedHash map[uint64][]*SeedBed, seedLen int, stepSize int, incomingFastqs <-chan *fastq.Fastq, outgoingSams chan<- *sam.SamAln, wg *sync.WaitGroup) {
	m, trace := swMatrixSetup(10000)
	for read := range incomingFastqs {
		//outgoingSams <- goGraphSmithWatermanMap(gg, read, seedHash, seedLen, stepSize, m, trace)
		outgoingSams <- GraphSmithWaterman(gg, read, seedHash, seedLen, stepSize, m, trace)
	}
	wg.Done()
}

func gswPairEnd(gg *SimpleGraph, seedHash map[uint64][]*SeedBed, seedLen int, stepSize int, incomingFastqs <-chan *fastq.PairedEnd, outgoingSams chan<- *sam.PairedSamAln, wg *sync.WaitGroup) {
	m, trace := swMatrixSetup(10000)
	for read := range incomingFastqs {
		outgoingSams <- PairedEndAlign(gg, read, seedHash, seedLen, stepSize, m, trace)
	}
	wg.Done()
}
