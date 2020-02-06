package simpleGraph

import (
	"github.com/vertgenlab/gonomics/fastq"
	"github.com/vertgenlab/gonomics/sam"
	"sync"
)

func gswWorker(gg *SimpleGraph, seedHash map[uint64][]*SeedBed, seedLen int, stepSize int, incomingFastqs <-chan *fastq.Fastq, outgoingSams chan<- *sam.SamAln, wg *sync.WaitGroup) {
	m, trace := swMatrixSetup(10000)

	for read := range incomingFastqs {
		outgoingSams <- goGraphSmithWatermanMap(gg, read, seedHash, seedLen, stepSize, m, trace)
	}

	wg.Done()
}

func PairedEndAlign(gg *SimpleGraph, readPair *fastq.PairedEnd, seedHash map[uint64][]*SeedBed, seedLen int, stepSize int) *sam.PairedSamAln {
	//m, trace := swMatrixSetup(10000)
	var mappedPair sam.PairedSamAln = sam.PairedSamAln{FwdSam: nil, RevSam: nil}
	mappedPair.FwdSam = GraphSmithWaterman(gg, readPair.Fwd, seedHash, seedLen, stepSize)
	mappedPair.RevSam = GraphSmithWaterman(gg, readPair.Rev, seedHash, seedLen, stepSize)
	return &mappedPair
}
