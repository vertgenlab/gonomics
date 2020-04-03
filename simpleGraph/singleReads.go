package simpleGraph

import (
	"github.com/vertgenlab/gonomics/fastq"
	"github.com/vertgenlab/gonomics/sam"
	"log"
	"sync"
)

func GswSingleReadWrap(ref *SimpleGraph, readOne string, output string, threads int, seedLen int, stepSize int, header *sam.SamHeader) {
	log.Printf("GSW!\n")
	log.Printf("Single end reads detected...\n")
	log.Printf("Aligning with the following settings: threads=%d, seedLen=%d, stepSize=%d\n\n", threads, seedLen, stepSize)
	log.Printf("Indexing the genome...\n")
	seedHash := indexGenomeIntoMap(ref.Nodes, seedLen, stepSize)
	var wgAlign, wgWrite sync.WaitGroup
	var scoreMatrix = HumanChimpTwoScoreMatrix
	log.Printf("Setting up read and write channels...\n")
	fastqPipe := make(chan *fastq.FastqBig, 824)
	samPipe := make(chan *sam.SamAln, 824)
	go fastq.ReadBigToChan(readOne, fastqPipe)
	wgAlign.Add(threads)
	for i := 0; i < threads; i++ {
		go gswWorkerMemPool(ref, seedHash, seedLen, stepSize, scoreMatrix, fastqPipe, samPipe, &wgAlign)
	}
	wgWrite.Add(1)
	go sam.SamChanToFile(samPipe, output, header, &wgWrite)
	wgAlign.Wait()
	close(samPipe)
	log.Printf("Aligners finished and channel closed\n")
	wgWrite.Wait()
	log.Printf("Sam writer finished and we are all done\n")
}
