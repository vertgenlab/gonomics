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
	settings := &GraphSettings{
		ScoreMatrix:    scoreMatrix,
		GapPenalty:     -600,
		OpenGapPenalty: -150,
		TileSize:       seedLen,
		StepSize:       stepSize,
	}
	allocateMemory := NewMemoryAllocation(defaultMatrixSize)
	seedBuildHelper := NewSeedBuilder()
	scorekeeper := scoreKeeper{}
	for read := range inputChan {
		outputChan <- *GraphSmithWatermanToGiraf(gg, read, seedHash, settings, allocateMemory, scorekeeper, seedBuildHelper)
	}
	wg.Done()
}

func RoutineFqPairToGiraf(gg *GenomeGraph, seedHash map[uint64][]uint64, seedLen int, stepSize int, scoreMatrix [][]int64, input <-chan fastq.PairedEndBig, output chan<- giraf.GirafPair, wg *sync.WaitGroup) {
	settings := &GraphSettings{
		ScoreMatrix:    scoreMatrix,
		GapPenalty:     -600,
		OpenGapPenalty: -150,
		TileSize:       seedLen,
		StepSize:       stepSize,
	}
	allocateMemory := NewMemoryAllocation(defaultMatrixSize)
	seedBuildHelper := NewSeedBuilder()
	scorekeeper := scoreKeeper{}
	for read := range input {
		output <- WrapPairGiraf(gg, read, seedHash, settings, allocateMemory, scorekeeper, seedBuildHelper)
	}
	wg.Done()
}

func RoutineGirafToSamSingle(gg *GenomeGraph, seedHash map[uint64][]uint64, seedLen int, stepSize int, scoreMatrix [][]int64, inputChan <-chan fastq.FastqBig, outputChan chan<- sam.Sam, wg *sync.WaitGroup) {
	settings := &GraphSettings{
		ScoreMatrix:    scoreMatrix,
		GapPenalty:     -600,
		OpenGapPenalty: -150,
		TileSize:       seedLen,
		StepSize:       stepSize,
	}

	allocateMemory := NewMemoryAllocation(defaultMatrixSize)

	seedBuildHelper := NewSeedBuilder()
	scorekeeper := scoreKeeper{
		leftAlignment:  make([]cigar.Cigar, 0, 1),
		rightAlignment: make([]cigar.Cigar, 0, 1),
	}
	for read := range inputChan {
		outputChan <- GirafToSam(GraphSmithWatermanToGiraf(gg, read, seedHash, settings, allocateMemory, scorekeeper, seedBuildHelper))
	}
	wg.Done()
}

func RoutineGirafToSam(gg *GenomeGraph, seedHash map[uint64][]uint64, seedLen int, stepSize int, scoreMatrix [][]int64, input <-chan fastq.PairedEndBig, output chan<- sam.Sam, wg *sync.WaitGroup) {
	settings := &GraphSettings{
		ScoreMatrix:    scoreMatrix,
		GapPenalty:     -600,
		OpenGapPenalty: -150,
		TileSize:       seedLen,
		StepSize:       stepSize,
	}
	allocateMemory := NewMemoryAllocation(defaultMatrixSize)
	scorekeeper := scoreKeeper{}
	seedBuildHelper := NewSeedBuilder()
	var pair sam.MatePair
	for read := range input {
		pair = GirafPairToSam(WrapPairGiraf(gg, read, seedHash, settings, allocateMemory, scorekeeper, seedBuildHelper))
		output <- pair.Fwd
		output <- pair.Rev
	}
	wg.Done()
}
