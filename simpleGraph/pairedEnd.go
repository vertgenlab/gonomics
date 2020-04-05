package simpleGraph

import (
	"fmt"
	"github.com/vertgenlab/gonomics/fastq"
	"github.com/vertgenlab/gonomics/sam"
	"log"
	"math"
	"sync"
	"time"
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
	log.Printf("Paired end reads detected...\n")

	log.Printf("Indexing the genome...\n\n")
	seedHash := indexGenomeIntoMap(ref.Nodes, seedLen, stepSize)
	var wgAlign, wgWrite sync.WaitGroup
	//log.Printf("Setting up read and write channels...\n\n")
	fastqPipe := make(chan *fastq.PairedEndBig, 824)
	samPipe := make(chan *sam.PairedSamAln, 824)
	go fastq.ReadPairBigToChan(readOne, readTwo, fastqPipe)
	wgWrite.Add(1)
	go sam.SamChanPairToFile(samPipe, output, header, &wgWrite)

	log.Printf("Scoring matrix used:\n%s\n", viewMatrix(scoreMatrix))
	log.Printf("Aligning with the following settings:\n\t\tthreads=%d, seedLen=%d, stepSize=%d\n\n", threads, seedLen, stepSize)
	wgAlign.Add(threads)
	time.Sleep(5 * time.Second)
	log.Printf("Aligning sequence to genome graph...")
	start := time.Now()
	for i := 0; i < threads; i++ {
		go gswWorkerPairedEnd(ref, seedHash, seedLen, stepSize, scoreMatrix, fastqPipe, samPipe, &wgAlign)
	}
	wgAlign.Wait()
	stop := time.Now()
	close(samPipe)
	wgWrite.Wait()
	log.Printf("GSW aligner finished in %.1f seconds\n", stop.Sub(start).Seconds())
	log.Printf("Enjoy analyzing your data!\n\n--xoxo GG\n")
}

func viewMatrix(m [][]int64) string {
	var message string = ""
	message += fmt.Sprintf("\t\t %d\t%d\t%d\t%d\n\t\t%d\t%d\t%d\t%d\n\t\t%d\t%d\t%d\t%d\n\t\t%d\t%d\t%d\t %d\n", m[0][0], m[0][1], m[0][2], m[0][3], m[1][0], m[1][1], m[1][2], m[1][3], m[2][0], m[2][1], m[2][2], m[2][3], m[3][0], m[3][1], m[3][2], m[3][3])
	return message
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
