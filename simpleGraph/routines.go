package simpleGraph

import (
	"github.com/vertgenlab/gonomics/fastq"
	"github.com/vertgenlab/gonomics/giraf"
	"github.com/vertgenlab/gonomics/sam"
	"sync"
)

//Goroutine worker functions
func RoutineFqToGiraf(gg *SimpleGraph, seedHash map[uint64][]uint64, seedLen int, stepSize int, scoreMatrix [][]int64, inputChan <-chan fastq.FastqBig, outputChan chan<- giraf.Giraf, wg *sync.WaitGroup) {
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

func RoutineFqPairToGiraf(gg *SimpleGraph, seedHash map[uint64][]uint64, seedLen int, stepSize int, scoreMatrix [][]int64, input <-chan fastq.PairedEndBig, output chan<- giraf.GirafPair, wg *sync.WaitGroup) {
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

func RoutineGirafToSamSingle(gg *SimpleGraph, seedHash map[uint64][]uint64, seedLen int, stepSize int, scoreMatrix [][]int64, inputChan <-chan fastq.FastqBig, outputChan chan<- *sam.SamAln, wg *sync.WaitGroup) {
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

func RoutineGirafToSam(gg *SimpleGraph, seedHash map[uint64][]uint64, seedLen int, stepSize int, scoreMatrix [][]int64, input <-chan fastq.PairedEndBig, output chan<- *sam.PairedSamAln, wg *sync.WaitGroup) {
	matrix := NewSwMatrix(defaultMatrixSize)
	seedPool := NewMemSeedPool()
	dnaPool := NewDnaPool()
	scorekeeper := scoreKeeper{}
	dynamicKeeper := dynamicScoreKeeper{}
	seedBuildHelper := newSeedBuilder()
	for read := range input {
		output <- GirafPairToSam(WrapPairGiraf(gg, read, seedHash, seedLen, stepSize, &matrix, scoreMatrix, &seedPool, &dnaPool, scorekeeper, dynamicKeeper, seedBuildHelper))
	}
	wg.Done()
}