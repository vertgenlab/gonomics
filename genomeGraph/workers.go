package genomeGraph

import (
	"sync"

	"github.com/vertgenlab/gonomics/cigar"
	"github.com/vertgenlab/gonomics/fastq"
	"github.com/vertgenlab/gonomics/sam"
)

// this is just for speed testing to see how much of the speed slowdown is due to alignment time.
func gswNothingWorker(gg *GenomeGraph, seedHash map[uint64][]uint64, seedLen int, stepSize int, incomingFastqs <-chan fastq.FastqBig, outgoingSams chan<- *sam.Sam, wg *sync.WaitGroup) {
	var curr *sam.Sam
	for read := range incomingFastqs {
		curr = &sam.Sam{QName: read.Name, Flag: 4, RName: "*", Pos: 0, MapQ: 255, Cigar: []cigar.Cigar{{Op: '*'}}, RNext: "*", PNext: 0, TLen: 0, Seq: read.Seq, Qual: string(read.Qual), Extra: "BZ:i:0"}
		outgoingSams <- curr
	}
	wg.Done()
}

// func gswWorkerMemPool(genome *GenomeGraph) {
// 	var tileSize int = 32
// 	var stepSize int = 32

// 	index := IndexGenomeIntoMap(genome.Nodes, tileSize, stepSize)

// 	var wg sync.WaitGroup
// 	var numWorkers int = 6
// 	input := make(chan fastq.FastqBig)

// 	matrix := NewSwMatrix(defaultMatrixSize)
// 	scores := HumanChimpTwoScoreMatrix

// 	seedPool := NewMemSeedPool()
// 	dnaPool := NewDnaPool()
// 	seedBuildHelper := newSeedBuilder()

// 	scorekeeper := scoreKeeper{}
// 	dynamicKeeper := dynamicScoreKeeper{}

// 	fastq.WritePair("testdata/simReads_R1.fq", "testdata/simReads_R2.fq")
// 	output := make(chan giraf.Giraf)

// 	go fastq.ReadPairBigToChan("testdata/simReads_R1.fq", "testdata/simReads_R2.fq", input)

// 	wg.Add(numWorkers)
// 	for i := 0; i < numWorkers; i++ {
// 		go GraphSmithWatermanToGiraf(genome, read, index, tileSize, stepSize, &matrix, scores, &seedPool, &dnaPool, scorekeeper, dynamicKeeper, seedBuildHelper)
// 	}
// 	go giraf.GirafPairChanToFile(output, samPipe, &wg)
// 	wg.Add(1)
// 	wg.Wait()
// 	close(samPipe)
// 	log.Printf("Aligners finished and channel closed\n")
// 	wg.Wait()
// 	log.Printf("Sam writer finished and we are all done\n")
// 	//stop := time.Now()
// 	wg.Done()
// }
