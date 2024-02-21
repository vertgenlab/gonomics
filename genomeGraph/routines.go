package genomeGraph

import (
	"sync"

	"github.com/vertgenlab/gonomics/cigar"
	"github.com/vertgenlab/gonomics/fastq"
	"github.com/vertgenlab/gonomics/giraf"
	"github.com/vertgenlab/gonomics/sam"
)

// Goroutine worker functions.
func RoutineFqToGiraf(gg *GenomeGraph, seedHash map[uint64][]uint64, seedLen int, stepSize int, scoreMatrix [][]int64, inputChan <-chan fastq.FastqBig, outputChan chan<- giraf.Giraf, wg *sync.WaitGroup) {
	matrix := NewSwMatrix(defaultMatrixSize)
	seedPool := NewMemSeedPool()
	dnaPool := NewDnaPool()
	seedBuildHelper := newSeedBuilder()
	scorekeeper := scoreKeeper{}
	dynamicKeeper := dynamicScoreKeeper{}
	for read := range inputChan {
		outputChan <- *GraphSmithWatermanToGiraf(gg, read, seedHash, seedLen, stepSize, &matrix, scoreMatrix, &seedPool, &dnaPool, scorekeeper, dynamicKeeper, seedBuildHelper)
	}
	wg.Done()
}

func RoutineFqPairToGiraf(gg *GenomeGraph, seedHash map[uint64][]uint64, seedLen int, stepSize int, scoreMatrix [][]int64, input <-chan fastq.PairedEndBig, output chan<- giraf.GirafPair, wg *sync.WaitGroup) {
	matrix := NewSwMatrix(defaultMatrixSize)
	seedPool := NewMemSeedPool()
	dnaPool := NewDnaPool()
	seedBuildHelper := newSeedBuilder()
	scorekeeper := scoreKeeper{}
	dynamicKeeper := dynamicScoreKeeper{}
	for read := range input {
		output <- WrapPairGiraf(gg, read, seedHash, seedLen, stepSize, &matrix, scoreMatrix, &seedPool, &dnaPool, scorekeeper, dynamicKeeper, seedBuildHelper)
	}
	wg.Done()
}

func RoutineGirafToSamSingle(gg *GenomeGraph, seedHash map[uint64][]uint64, seedLen int, stepSize int, scoreMatrix [][]int64, inputChan <-chan fastq.FastqBig, outputChan chan<- sam.Sam, wg *sync.WaitGroup) {
	matrix := NewSwMatrix(defaultMatrixSize)
	seedPool := NewMemSeedPool()
	dnaPool := NewDnaPool()
	seedBuildHelper := newSeedBuilder()
	scorekeeper := scoreKeeper{}
	dynamicKeeper := dynamicScoreKeeper{}
	for read := range inputChan {
		outputChan <- GirafToSam(GraphSmithWatermanToGiraf(gg, read, seedHash, seedLen, stepSize, &matrix, scoreMatrix, &seedPool, &dnaPool, scorekeeper, dynamicKeeper, seedBuildHelper))
	}
	wg.Done()
}

func RoutineGirafToSam(gg *GenomeGraph, seedHash map[uint64][]uint64, seedLen int, stepSize int, scoreMatrix [][]int64, input <-chan fastq.PairedEndBig, output chan<- sam.Sam, wg *sync.WaitGroup) {
	matrix := NewSwMatrix(defaultMatrixSize)
	seedPool := NewMemSeedPool()
	dnaPool := NewDnaPool()
	scorekeeper := scoreKeeper{}
	dynamicKeeper := dynamicScoreKeeper{}
	seedBuildHelper := newSeedBuilder()
	var pair sam.MatePair
	for read := range input {
		pair = GirafPairToSam(WrapPairGiraf(gg, read, seedHash, seedLen, stepSize, &matrix, scoreMatrix, &seedPool, &dnaPool, scorekeeper, dynamicKeeper, seedBuildHelper))
		output <- pair.Fwd
		output <- pair.Rev
	}
	wg.Done()
}

// this is just for speed testing to see how much of the speed slowdown is due to alignment time.
func gswNothingWorker(gg *GenomeGraph, seedHash map[uint64][]uint64, seedLen int, stepSize int, incomingFastqs <-chan fastq.FastqBig, outgoingSams chan<- *sam.Sam, wg *sync.WaitGroup) {
	var curr *sam.Sam
	for read := range incomingFastqs {
		curr = &sam.Sam{QName: read.Name, Flag: 4, RName: "*", Pos: 0, MapQ: 255, Cigar: []cigar.Cigar{{Op: '*'}}, RNext: "*", PNext: 0, TLen: 0, Seq: read.Seq, Qual: string(read.Qual), Extra: "BZ:i:0"}
		outgoingSams <- curr
	}
	wg.Done()
}
