package simpleGraph

import (
	"github.com/vertgenlab/gonomics/fastq"
	"github.com/vertgenlab/gonomics/sam"
	"log"
	"math"
	"sync"
)

func PairedEndTwoBitAlign(gg *SimpleGraph, readPair *fastq.PairedEndBig, seedHash map[uint64][]uint64, seedLen int, stepSize int, scoreMatrix [][]int64, m [][]int64, trace [][]rune, memoryPool **SeedDev) *sam.PairedSamAln {
	var mappedPair sam.PairedSamAln = sam.PairedSamAln{FwdSam: nil, RevSam: nil}
	mappedPair.FwdSam = GraphSmithWaterman(gg, readPair.Fwd, seedHash, seedLen, stepSize, scoreMatrix, m, trace, memoryPool)
	mappedPair.RevSam = GraphSmithWaterman(gg, readPair.Rev, seedHash, seedLen, stepSize, scoreMatrix, m, trace, memoryPool)
	mappedPair.FwdSam.Flag += 64
	mappedPair.RevSam.Flag += 128
	if math.Abs(float64(mappedPair.FwdSam.Pos-mappedPair.RevSam.Pos)) < 10000 {
		mappedPair.FwdSam.Flag += 2
		mappedPair.RevSam.Flag += 2
	}
	return &mappedPair
}

func gswWorkerPairedEnd(gg *SimpleGraph, seedHash map[uint64][]uint64, seedLen int, stepSize int, scoreMatrix [][]int64, incomingFastqs <-chan *fastq.PairedEndBig, outgoingSams chan<- *sam.PairedSamAln, wg *sync.WaitGroup) {
	m, trace := swMatrixSetup(10000)
	memChunk := make([]SeedDev, 100000)
	for i := 0; i < len(memChunk)-1; i++ {
		memChunk[i].Next = &memChunk[i+1]
	}
	memStart := &(memChunk[0])
	for read := range incomingFastqs {
		outgoingSams <- PairedEndTwoBitAlign(gg, read, seedHash, seedLen, stepSize, scoreMatrix, m, trace, &memStart)
	}
	wg.Done()
}

func GSWsPair(ref *SimpleGraph, readOne string, readTwo string, output string, threads int, seedLen int, stepSize int, scoreMatrix [][]int64, header *sam.SamHeader) {
	log.SetFlags(log.Ldate | log.Ltime)
	log.Printf("GSW!\n")
	log.Printf("Paired end reads detected...\n")
	log.Printf("Aligning with the following settings: threads=%d, seedLen=%d, stepSize=%d\n\n", threads, seedLen, stepSize)
	log.Printf("Indexing the genome...\n")
	seedHash := indexGenomeIntoMap(ref.Nodes, seedLen, stepSize)
	var wgAlign, wgWrite sync.WaitGroup
	log.Printf("Setting up read and write channels...\n")
	fastqPipe := make(chan *fastq.PairedEndBig, 824)
	samPipe := make(chan *sam.PairedSamAln, 824)
	go fastq.ReadPairBigToChan(readOne, readTwo, fastqPipe)
	wgAlign.Add(threads)
	for i := 0; i < threads; i++ {
		go gswWorkerPairedEnd(ref, seedHash, seedLen, stepSize, scoreMatrix, fastqPipe, samPipe, &wgAlign)
	}
	wgWrite.Add(1)
	go sam.SamChanPairToFile(samPipe, output, header, &wgWrite)
	wgAlign.Wait()
	close(samPipe)
	log.Printf("Aligners finished and channel closed\n")
	wgWrite.Wait()
	log.Printf("Sam writer finished and we are all done\n")
}

func PairedEndAlign(gg *SimpleGraph, readPair *fastq.PairedEnd, seedHash map[uint64][]*SeedBed, seedLen int, stepSize int, m [][]int64, trace [][]rune) *sam.PairedSamAln {
	var mappedPair sam.PairedSamAln = sam.PairedSamAln{FwdSam: nil, RevSam: nil}
	mappedPair.FwdSam = oldGraphSmithWaterman(gg, readPair.Fwd, seedHash, seedLen, stepSize, m, trace)
	mappedPair.RevSam = oldGraphSmithWaterman(gg, readPair.Rev, seedHash, seedLen, stepSize, m, trace)
	mappedPair.FwdSam.Flag += 64
	mappedPair.RevSam.Flag += 128
	if math.Abs(float64(mappedPair.FwdSam.Pos-mappedPair.RevSam.Pos)) < 10000 {
		mappedPair.FwdSam.Flag += 2
		mappedPair.RevSam.Flag += 2
	}
	return &mappedPair
}
