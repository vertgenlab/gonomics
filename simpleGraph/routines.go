package simpleGraph

import (
	"github.com/vertgenlab/gonomics/fastq"
	"github.com/vertgenlab/gonomics/sam"
)

func gswWorker(gg *SimpleGraph, seedHash map[uint64][]*SeedBed, seedLen int, stepSize int, incomingFastqs <-chan *fastq.Fastq, outgoingSams chan<- *sam.SamAln) {
	m, trace := swMatrixSetup(10000)

	for read := range incomingFastqs {
		outgoingSams <- goGraphSmithWatermanMap(gg, read, seedHash, seedLen, stepSize, m, trace)
	}
}
