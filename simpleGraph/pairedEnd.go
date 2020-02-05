package simpleGraph

import (
	"github.com/vertgenlab/gonomics/fastq"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/sam"
	"log"
	"os"
	"sync"
)

func GSWsBatchPair(ref *SimpleGraph, readOne string, readTwo string, output string, groupSize int) {

	var seedLen int = 30
	var stepSize int = 29
	log.Printf("Indexing the genome...\n")
	seedHash := IndexGenomeIntoMap(ref.Nodes, seedLen, stepSize)
	var fq *fastq.PairedEnd
	var groups []*fastq.PairedEnd = make([]*fastq.PairedEnd, 0, groupSize)
	var done bool
	//fastq input
	file := fileio.EasyOpen(readOne)
	defer file.Close()

	file2 := fileio.EasyOpen(readTwo)
	defer file.Close()

	samRecords, _ := os.Create(output)
	defer samRecords.Close()

	header := NodesHeader(ref.Nodes)
	sam.WriteHeaderToFileHandle(samRecords, header)
	batchSize := 5
	batches := make([][]*fastq.PairedEnd, 0, batchSize)
	var batchID int = 0
	var wg sync.WaitGroup
	for fq, done = fastq.NextFastqPair(file, file2); !done; fq, done = fastq.NextFastqPair(file, file2) {
		groups = append(groups, fq)

		if len(groups) == groupSize {
			//send off group

			batches = append(batches, groups)
			groups = make([]*fastq.PairedEnd, 0, groupSize)

			if len(batches) == batchSize {
				log.Printf("Aligning reads...\n")
				alignBatchPairGroup(ref, batches, seedHash, seedLen, stepSize, 824, samRecords, batchID, &wg)
				batches = make([][]*fastq.PairedEnd, 0, batchSize)
			}
			//zero out
			//alignGroups(ref, groups, seedHash, seedLen, 8024, samRecords)

			//go routineGenomeGraph(ref, <-batches, seedHash, seedLen, 8204, &wg, samRecords, batchID)

		}
	}
	//wg.Wait()
	//last case check length groups
	if len(batches) != 0 {
		if len(groups) != 0 {
			batches = append(batches, groups)
		}

		alignBatchPairGroup(ref, batches, seedHash, seedLen, stepSize, 824, samRecords, batchID, &wg)
	}
	//wg.Wait()
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
			//m, trace := SwMatrixSetup(10000)
			for fqs := range input {
				//mappedPair = GraphSmithWaterman(gg, fqs.Fwd, seedHash, seedLen, stepSize, m, trace)
				//mappedPair = GraphSmithWaterman(gg, fqs.Rev, seedHash, seedLen, stepSize, m, trace)
				output <- PairedEndAlign(gg, fqs, seedHash, seedLen, stepSize)
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
