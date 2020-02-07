package simpleGraph

import (
	"github.com/vertgenlab/gonomics/fastq"
	//"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/sam"
	"log"
	"os"
	"sync"
)

func PairedEndAlign(gg *SimpleGraph, readPair *fastq.PairedEnd, seedHash map[uint64][]*SeedBed, seedLen int, stepSize int, m [][]int64, trace [][]rune) *sam.PairedSamAln {
	//m, trace := swMatrixSetup(10000)
	var mappedPair sam.PairedSamAln = sam.PairedSamAln{FwdSam: nil, RevSam: nil}
	mappedPair.FwdSam = GraphSmithWaterman(gg, readPair.Fwd, seedHash, seedLen, stepSize, m, trace)
	mappedPair.RevSam = GraphSmithWaterman(gg, readPair.Rev, seedHash, seedLen, stepSize, m, trace)
	return &mappedPair
}
//TODO: Does not work on graph with edges. Work in progress
func GSWsBatchPair(ref *SimpleGraph, readOne string, readTwo string, output string) {

	var seedLen int = 32
	var stepSize int = seedLen - 1
	var numWorkers int = 8
	log.Printf("Indexing the genome...\n")
	seedHash := IndexGenomeIntoMap(ref.Nodes, seedLen, stepSize)

	samRecords, _ := os.Create(output)
	defer samRecords.Close()
	header := NodesHeader(ref.Nodes)
	sam.WriteHeaderToFileHandle(samRecords, header)

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
	wgAlign.Add(1)

	go sam.SamChanPairToFile(samPipe, samRecords, &wgWrite)
	wgAlign.Wait()
	close(samPipe)
	log.Printf("Aligners finished and channel closed\n")
	wgWrite.Wait()
	log.Printf("Sam writer finished and we are all done\n")
}

func routineGswReadPair(gg *SimpleGraph, reads []*fastq.PairedEnd, seedHash map[uint64][]*SeedBed, seedLen int, stepSize int, workers int, wg *sync.WaitGroup, file *os.File, batchNum int) {

	batchId := batchNum
	numJobs := len(reads)
	input := make(chan *fastq.PairedEnd, numJobs)
	output := make(chan *sam.PairedSamAln, numJobs)
	//var mappedPair = make([]*sam.SamAln, 2)
	for w := 0; w < workers; w++ {
		//wg.Add(1)
		go func(gg *SimpleGraph, reads []*fastq.PairedEnd, seedHash map[uint64][]*SeedBed, seedLen int, input <-chan *fastq.PairedEnd, output chan<- *sam.PairedSamAln) {
			//map reads using worker pool
			m, trace := SwMatrixSetup(10000)
			for fqs := range input {
				//mappedPair = GraphSmithWaterman(gg, fqs.Fwd, seedHash, seedLen, stepSize, m, trace)
				//mappedPair = GraphSmithWaterman(gg, fqs.Rev, seedHash, seedLen, stepSize, m, trace)
				output <- PairedEndAlign(gg, fqs, seedHash, seedLen, stepSize, m, trace)
			}

		}(gg, reads, seedHash, seedLen, input, output)
	}

	for i := 0; i < numJobs; i++ {
		//send reads off to be aligned
		input <- reads[i]
	}
	close(input)
	for j := 0; j < numJobs; j++ {
		//print something here
		//or write to file
		sam.WriteAlnPairToFileHandle(file, <-output)
		log.Printf("Finished aligning %d reads in batch %d\n", j, batchId)
		//log.Printf("%s\n", sam.SamAlnToString(<-output))

	}
	wg.Done()
}

func alignBatchPairGroup(gg *SimpleGraph, batch [][]*fastq.PairedEnd, seedHash map[uint64][]*SeedBed, seedLen int, stepSize int, workers int, file *os.File, batchNum int, wg *sync.WaitGroup) {
	for i := 0; i < len(batch); i++ {
		wg.Add(1)
		batchNum++
		go routineGswReadPair(gg, batch[i], seedHash, seedLen, stepSize, 824, wg, file, batchNum)
	}
	wg.Wait()
}
