package genomeGraph

import (
	"github.com/vertgenlab/gonomics/cigar"
	"github.com/vertgenlab/gonomics/fastq"
	"github.com/vertgenlab/gonomics/sam"
	"sync"
)

// this is just for speed testing to see how much of the speed slowdown is due to alignment time
func gswNothingWorker(gg *GenomeGraph, seedHash map[uint64][]uint64, seedLen int, stepSize int, incomingFastqs <-chan *fastq.FastqBig, outgoingSams chan<- *sam.SamAln, wg *sync.WaitGroup) {
	var curr *sam.SamAln
	for read := range incomingFastqs {
		curr = &sam.SamAln{QName: read.Name, Flag: 4, RName: "*", Pos: 0, MapQ: 255, Cigar: []*cigar.Cigar{&cigar.Cigar{Op: '*'}}, RNext: "*", PNext: 0, TLen: 0, Seq: read.Seq, Qual: string(read.Qual), Extra: "BZ:i:0"}
		outgoingSams <- curr
	}
	wg.Done()
}

func gswWorkerMemPool(gg *GenomeGraph, seedHash map[uint64][]uint64, seedLen int, stepSize int, scoreMatrix [][]int64, incomingFastqs <-chan *fastq.FastqBig, outgoingSams chan<- *sam.SamAln, wg *sync.WaitGroup) {
	m, trace := swMatrixSetup(10000)
	for read := range incomingFastqs {
		outgoingSams <- GraphSmithWatermanMemPool(gg, read, seedHash, seedLen, stepSize, scoreMatrix, m, trace, nil)
	}
	wg.Done()
}
