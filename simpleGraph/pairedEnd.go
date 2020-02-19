package simpleGraph

import (
	"github.com/vertgenlab/gonomics/fastq"
	"github.com/vertgenlab/gonomics/sam"
	"log"
	"sync"
)

func PairedEndAlign(gg *SimpleGraph, readPair *fastq.PairedEnd, seedHash map[uint64][]*SeedBed, seedLen int, stepSize int, m [][]int64, trace [][]rune, noNs bool) *sam.PairedSamAln {
	//m, trace := swMatrixSetup(10000)
	var mappedPair sam.PairedSamAln = sam.PairedSamAln{FwdSam: nil, RevSam: nil}
	mappedPair.FwdSam = GraphSmithWaterman(gg, readPair.Fwd, seedHash, seedLen, stepSize, m, trace, noNs)
	mappedPair.RevSam = GraphSmithWaterman(gg, readPair.Rev, seedHash, seedLen, stepSize, m, trace, noNs)
	return &mappedPair
}

//TODO: Does not work on graph with edges. Work in progress
func GSWsBatchPair(ref *SimpleGraph, seedLen int, stepSize int, noNs bool, readOne string, readTwo string, output string) {

	//var seedLen int = 32
	//var stepSize int = seedLen - 1
	var numWorkers int = 8
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
		go gswPairEnd(ref, seedHash, seedLen, stepSize, fastqPipe, samPipe, &wgAlign, noNs)
	}
	wgAlign.Add(1)
	go sam.SamChanPairToFile(samPipe, output, header, &wgWrite)

	wgAlign.Wait()
	close(samPipe)
	log.Printf("Aligners finished and channel closed\n")
	wgWrite.Wait()
	log.Printf("Sam writer finished and we are all done\n")
}
