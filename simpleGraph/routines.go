package simpleGraph

import (
	"github.com/vertgenlab/gonomics/fastq"
	"github.com/vertgenlab/gonomics/giraf"
	"github.com/vertgenlab/gonomics/sam"
	"sync"
)

//Goroutine worker functions
func RoutineFqToGiraf(gg *SimpleGraph, seedHash map[uint64][]uint64, seedLen int, stepSize int, scoreMatrix [][]int64, inputChan <-chan *fastq.FastqBig, outputChan chan<- *giraf.Giraf, wg *sync.WaitGroup) {
	m, trace := swMatrixSetup(10000)
	/*memChunk := make([]SeedDev, 100000)
	for i := 0; i < len(memChunk)-1; i++ {
		memChunk[i].Next = &memChunk[i+1]
	}
	memStart := &(memChunk[0])*/
	for read := range inputChan {
		outputChan <- GraphSmithWatermanToGiraf(gg, read, seedHash, seedLen, stepSize, scoreMatrix, m, trace)
	}
	wg.Done()
}

func RoutineFqPairToGiraf(gg *SimpleGraph, seedHash map[uint64][]uint64, seedLen int, stepSize int, scoreMatrix [][]int64, input <-chan *fastq.PairedEndBig, output chan<- *giraf.GirafPair, wg *sync.WaitGroup) {
	m, trace := swMatrixSetup(10000)
	/*memChunk := make([]SeedDev, 100000)
	for i := 0; i < len(memChunk)-1; i++ {
		memChunk[i].Next = &memChunk[i+1]
	}
	memStart := &(memChunk[0])*/
	for read := range input {
		output <- WrapPairGiraf(gg, read, seedHash, seedLen, stepSize, scoreMatrix, m, trace)
	}
	wg.Done()
}

func RoutineGirafToSamSingle(gg *SimpleGraph, seedHash map[uint64][]uint64, seedLen int, stepSize int, scoreMatrix [][]int64, inputChan <-chan *fastq.FastqBig, outputChan chan<- *sam.SamAln, wg *sync.WaitGroup) {
	m, trace := swMatrixSetup(10000)
	/*memChunk := make([]SeedDev, 100000)
	for i := 0; i < len(memChunk)-1; i++ {
		memChunk[i].Next = &memChunk[i+1]
	}
	memStart := &(memChunk[0])*/
	for read := range inputChan {
		outputChan <- GirafToSam(GraphSmithWatermanToGiraf(gg, read, seedHash, seedLen, stepSize, scoreMatrix, m, trace))
	}
	wg.Done()
}

func RoutineGirafToSam(gg *SimpleGraph, seedHash map[uint64][]uint64, seedLen int, stepSize int, scoreMatrix [][]int64, inputChan <-chan *fastq.PairedEndBig, outputChan chan<- *sam.PairedSamAln, wg *sync.WaitGroup) {
	m, trace := swMatrixSetup(10000)
	/*memChunk := make([]SeedDev, 100000)
	for i := 0; i < len(memChunk)-1; i++ {
		memChunk[i].Next = &memChunk[i+1]
	}
	memStart := &(memChunk[0])*/
	for read := range inputChan {
		outputChan <- GirafLiftoverToSam(gg, read, seedHash, seedLen, stepSize, scoreMatrix, m, trace)
	}
	wg.Done()
}
func gswWorker(gg *SimpleGraph, seedHash map[uint64][]*SeedBed, seedLen int, stepSize int, incomingFastqs <-chan *fastq.Fastq, outgoingSams chan<- *sam.SamAln, wg *sync.WaitGroup) {
	m, trace := swMatrixSetup(10000)
	for read := range incomingFastqs {
		//outgoingSams <- devGraphSmithWaterman(gg, read, seedHash, seedLen, stepSize, m, trace, HumanChimpTwoScoreMatrix)
		outgoingSams <- oldGraphSmithWaterman(gg, read, seedHash, seedLen, stepSize, m, trace)
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
