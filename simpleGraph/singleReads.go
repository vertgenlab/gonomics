package simpleGraph

import (
	"github.com/vertgenlab/gonomics/chromInfo"
	"github.com/vertgenlab/gonomics/fastq"
	"github.com/vertgenlab/gonomics/sam"
	"log"
	"sync"
)

func GswSingleReadWrap(ref *SimpleGraph, readOne string, output string, threads int, seedLen int, stepSize int, chrSize map[string]*chromInfo.ChromInfo) {
	//var seedLen int = kMer
	//var stepSize int = seedLen - 1
	log.Printf("Reading reference...\n")
	//ref, chrSize := Read(filename)
	var numWorkers int = threads
	log.Printf("Indexing the genome...\n")
	seedHash := IndexGenomeIntoMap(ref.Nodes, seedLen, stepSize)
	header := sam.ChromInfoMapSamHeader(chrSize)
	var wgAlign, wgWrite sync.WaitGroup

	log.Printf("Setting up goroutine channels...\n")
	fastqPipe := make(chan *fastq.Fastq, 824)

	//log.Printf("Making sam channel...\n")
	samPipe := make(chan *sam.SamAln, 824)
	go fastq.ReadToChan(readOne, fastqPipe)
	wgAlign.Add(numWorkers)
	log.Printf("Aligning fastqs to graph...\n")
	for i := 0; i < numWorkers; i++ {
		go gswWorker(ref, seedHash, seedLen, stepSize, fastqPipe, samPipe, &wgAlign)
	}
	wgWrite.Add(1)

	go sam.SamChanToFile(samPipe, output, header, &wgWrite)
	wgAlign.Wait()
	close(samPipe)
	wgWrite.Wait()
	log.Printf("Finished aligning fastqs!!\n")
}
