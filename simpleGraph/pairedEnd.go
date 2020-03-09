package simpleGraph

import (
	"github.com/vertgenlab/gonomics/fastq"
	"github.com/vertgenlab/gonomics/sam"
	"log"
	"math"
	"sync"
)

func PairedEndAlign(gg *SimpleGraph, readPair *fastq.PairedEnd, seedHash map[uint64][]*SeedBed, seedLen int, stepSize int, m [][]int64, trace [][]rune) *sam.PairedSamAln {
	//m, trace := swMatrixSetup(10000)
	var mappedPair sam.PairedSamAln = sam.PairedSamAln{FwdSam: nil, RevSam: nil}
	mappedPair.FwdSam = GraphSmithWaterman(gg, readPair.Fwd, seedHash, seedLen, stepSize, m, trace)
	mappedPair.RevSam = GraphSmithWaterman(gg, readPair.Rev, seedHash, seedLen, stepSize, m, trace)
	mappedPair.FwdSam.Flag += 64
	mappedPair.RevSam.Flag += 128
	if math.Abs(float64(mappedPair.FwdSam.Pos-mappedPair.RevSam.Pos)) < 10000 {
		mappedPair.FwdSam.Flag += 2
		mappedPair.RevSam.Flag += 2
	}
	return &mappedPair
}

func GSWsBatchPair(ref *SimpleGraph, readOne string, readTwo string, output string, threads int, seedLen int) {

	//var seedLen int = kMer
	var stepSize int = seedLen - 1
	var numWorkers int = threads
	log.Printf("Indexing the genome...\n")
	seedHash := IndexGenomeIntoMap(ref.Nodes, seedLen, stepSize)
	header := NodesHeader(ref.Nodes)
	var wgAlign, wgWrite sync.WaitGroup

	log.Printf("Making fastq channel...\n")
	fastqPipe := make(chan *fastq.PairedEnd, 824)

	log.Printf("Making sam channel...\n")
	samPipe := make(chan *sam.PairedSamAln, 824)
	go fastq.PairEndToChan(readOne, readTwo, fastqPipe)

	wgAlign.Add(numWorkers)
	for i := 0; i < numWorkers; i++ {
		go gswPairEnd(ref, seedHash, seedLen, stepSize, fastqPipe, samPipe, &wgAlign)
	}
	wgWrite.Add(1)
	go sam.SamChanPairToFile(samPipe, output, header, &wgWrite)
	wgAlign.Wait()
	close(samPipe)
	log.Printf("Aligners finished and channel closed\n")
	wgWrite.Wait()
	log.Printf("Sam writer finished and we are all done\n")
}

/*
func getTimeToAlign(startTime time.Time, done <-chan bool, numberAligned int, tick time.Ticker) {
	stopTime := time.Now()
	duration := stopTime.Sub(startTime)
	go func() {
		for {
			select {
			case <-done:
				return
			case <-tick.C:
				log.Printf("Aligned %d reads in %d\n", numberAligned, duration)
			}
		}
	}()
}

func WriteWithTime(incomingSams <-chan *sam.PairedSamAln, filename string, header *sam.SamHeader, wg *sync.WaitGroup) {
	file, _ := os.Create(filename)
	defer file.Close()
	sam.WriteHeaderToFileHandle(file, header)
	ticker := time.NewTicker(2 * time.Minute)
	done := make(chan bool)
	start := time.Now()
	var numFinished int = 0
	go getTimeToAlign(start, done, numFinished, *ticker)
	for alignedRead := range incomingSams {
		numFinished++
		sam.WriteAlnPairToFileHandle(file, alignedRead)
	}
	ticker.Stop()
	done <- true
	wg.Done()
}
*/
