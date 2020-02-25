package simpleGraph

/*
import (
	"github.com/vertgenlab/gonomics/fastq"
	//"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/sam"
	"log"
	"os"
	"sync"
	//"time"
)

func GSWsBatchDraft(ref *SimpleGraph, readOne string, output string) {

	var seedLen int = 32
	var stepSize int = seedLen - 1
	var numWorkers int = 16
	log.Printf("Indexing the genome...\n")
	seedHash := IndexGenomeIntoMap(ref.Nodes, seedLen, stepSize)

	samRecords, _ := os.Create(output)
	defer samRecords.Close()
	header := NodesHeader(ref.Nodes)
	sam.WriteHeaderToFileHandle(samRecords, header)

	var wgAlign, wgWrite sync.WaitGroup

	log.Printf("Making fastq channel...\n")
	fastqPipe := make(chan *fastq.Fastq, 824)

	log.Printf("Making sam channel...\n")
	samPipe := make(chan *sam.SamAln, 824)
	go fastq.ReadToChan(readOne, fastqPipe)

	wgAlign.Add(numWorkers)

	for i := 0; i < numWorkers; i++ {
		go gswWorker(ref, seedHash, seedLen, stepSize, fastqPipe, samPipe, &wgAlign)
	}
	wgAlign.Add(1)

	go sam.TestSamChanToFile(samPipe, samRecords, &wgWrite)
	wgAlign.Wait()
	close(samPipe)
	log.Printf("Aligners finished and channel closed\n")
	wgWrite.Wait()
	log.Printf("Sam writer finished and we are all done\n")
}*/

/*
func GSWsBatch(ref *SimpleGraph, input string, output string, groupSize int) {

	var seedLen int = 30
	var stepSize int = 24
	log.Printf("Indexing the genome...\n")
	seedHash := IndexGenomeIntoMap(ref.Nodes, seedLen, stepSize)
	var fq *fastq.Fastq
	var groups []*fastq.Fastq = make([]*fastq.Fastq, 0, groupSize)
	var done bool
	//fastq input
	file := fileio.EasyOpen(input)
	defer file.Close()

	samRecords, _ := os.Create(output)
	defer samRecords.Close()

	header := NodesHeader(ref.Nodes)
	sam.WriteHeaderToFileHandle(samRecords, header)

	for fq, done = fastq.NextFastq(file); !done; fq, done = fastq.NextFastq(file) {
		groups = append(groups, fq)
		if len(groups) == groupSize {
			//send off group
			alignGroups(ref, groups, seedHash, seedLen, stepSize, samRecords)
			//	groups = make([]*fastq.Fastq, 0, groupSize)
		}
	}
	if len(groups) != 0 {
		alignGroups(ref, groups, seedHash, seedLen, stepSize, samRecords)
	}
}

func alignGroups(gg *SimpleGraph, reads []*fastq.Fastq, seedHash map[uint64][]*SeedBed, seedLen int, stepSize int, file *os.File) {
	alignments := make([]*sam.SamAln, len(reads))
	log.Printf("Making fastq channel...\n")
	fastqPipe := make(chan *fastq.Fastq, len(reads))

	log.Printf("Making sam channel...\n")
	samPipe := make(chan *sam.SamAln, len(reads))

	//go gswWorker(gg, seedHash, seedLen, stepSize, fastqPipe, samPipe)
	//log.Printf("Waiting for 10 seconds and then aligning reads...\n")
	//time.Sleep(10 * time.Second)

	for i := 0; i < len(reads); i++ {
		fastqPipe <- reads[i]
	}
	for j := 0; j < len(reads); j++ {
		alignments[j] = <-samPipe
		//sam.WriteAlnToFileHandle(file, <-samPipe)
		//log.Printf("Finished aligning %d reads\n", j)
	}
	for records := range alignments {
		log.Printf("%s\n", sam.SamAlnToString(alignments[records]))
		//sam.WriteAlnToFileHandle(file, alignments[records])
	}
	//var mappedRead *sam.SamAln
	//for w := 0; w < workers; w++ {
	//wg.Add(1)
	//go func(gg *SimpleGraph, reads []*fastq.Fastq, seedHash map[uint64][]*SeedBed, seedLen int, input <-chan *fastq.Fastq, output chan<- *sam.SamAln) {
	//m, trace := swMatrixSetup(10000)
	//map reads using worker pool
	//for fq := range input {
	//	output <- GraphSmithWaterman(gg, fq, seedHash, seedLen, stepSize)
	//}

	//}(gg, reads, seedHash, seedLen, input, output)
	//}

	//for i := 0; i < numJobs; i++ {
	//send reads off to be aligned
	//	input <- reads[i]
	//}
	//close(input)
	//for j := 0; j < numJobs; j++ {
	//print something here
	//or write to file

	//log.Printf("%s\n", sam.SamAlnToString(<-output))

	//}
}

func routineGenomeGraph(gg *SimpleGraph, reads []*fastq.Fastq, seedHash map[uint64][]*SeedBed, seedLen int, stepSize int, workers int, wg *sync.WaitGroup, file *os.File, batchNum int) {

	batchId := batchNum
	numJobs := len(reads)
	input := make(chan *fastq.Fastq, numJobs)
	output := make(chan *sam.SamAln, numJobs)

	//go gswWorker(gg, seedHash, seedLen, stepSize, input, output)
	//var mappedRead *sam.SamAln
	//for w := 0; w < workers; w++ {
	//wg.Add(1)

	//	go func(gg *SimpleGraph, reads []*fastq.Fastq, seedHash map[uint64][]*SeedBed, seedLen int, input <-chan *fastq.Fastq, output chan<- *sam.SamAln) {
	//map reads using worker pool

	//		for fq := range input {
	//m, trace := swMatrixSetup(10000)
	//			output <- GraphSmithWaterman(gg, fq, seedHash, seedLen, stepSize)
	//		}

	//	}(gg, reads, seedHash, seedLen, input, output)
	//}

	for i := 0; i < numJobs; i++ {
		//send reads off to be aligned
		input <- reads[i]
	}
	close(input)
	for j := 0; j < numJobs; j++ {
		//print something here
		//or write to file
		sam.WriteAlnToFileHandle(file, <-output)
		log.Printf("Finished aligning %d reads in batch %d\n", j, batchId)
		//log.Printf("%s\n", sam.SamAlnToString(<-output))

	}
	wg.Done()
}

func alignBatchGroup(gg *SimpleGraph, batch [][]*fastq.Fastq, seedHash map[uint64][]*SeedBed, seedLen int, stepSize int, workers int, file *os.File, batchNum int, wg *sync.WaitGroup) {
	for i := 0; i < len(batch); i++ {
		wg.Add(1)
		batchNum++
		go routineGenomeGraph(gg, batch[i], seedHash, seedLen, stepSize, 824, wg, file, batchNum)
	}
	wg.Wait()
}

//GO Pher
//number of jobs is the length of reads
func TestGopher(gg *SimpleGraph, reads []*fastq.Fastq, seedHash map[uint64][]*SeedBed, seedLen int, stepSize int, workers int, wg *sync.WaitGroup) {
	numJobs := len(reads)
	input := make(chan *fastq.Fastq, numJobs)
	output := make(chan *sam.SamAln, numJobs)
	//var mappedRead *sam.SamAln
	for w := 0; w < workers; w++ {
		//wg.Add(1)
		go func(gg *SimpleGraph, reads []*fastq.Fastq, seedHash map[uint64][]*SeedBed, seedLen int, input <-chan *fastq.Fastq, output chan<- *sam.SamAln) {
			m, trace := swMatrixSetup(10000)
			//map reads using worker pool
			for fq := range input {
				output <- GraphSmithWaterman(gg, fq, seedHash, stepSize, seedLen, m, trace)
			}

		}(gg, reads, seedHash, seedLen, input, output)
	}

	for i := 0; i < numJobs; i++ {
		//send reads off to be aligned
		input <- reads[i]
	}
	close(input)
	for j := 0; j < numJobs; j++ {
		//<-output
		log.Printf("%s\n", sam.SamAlnToString(<-output))
		//print something here
		//or write to file
	}
	wg.Done()
}

func TestGroups(gg *SimpleGraph, reads []*fastq.Fastq, seedHash map[uint64][]*SeedBed, seedLen int, stepSize int, wg *sync.WaitGroup) {
	var mappedReads []*sam.SamAln = make([]*sam.SamAln, len(reads))
	for i := 0; i < len(reads); i++ {
		//wrapNoChan(gg, reads[i], seedHash, seedLen)
		m, trace := swMatrixSetup(10000)
		mappedReads[i] = GraphSmithWaterman(gg, reads[i], seedHash, stepSize, seedLen, m, trace)
	}
	for j := 0; j < len(mappedReads); j++ {
		log.Printf("%s\n", sam.SamAlnToString(mappedReads[j]))
	}
	wg.Done()
}

func GSWsBatches(ref *SimpleGraph, input string, output string, groupSize int) {

	var seedLen int = 30
	var stepSize int = 29
	log.Printf("Indexing the genome...\n")
	seedHash := IndexGenomeIntoMap(ref.Nodes, seedLen, stepSize)
	var fq *fastq.Fastq
	var groups []*fastq.Fastq = make([]*fastq.Fastq, 0, groupSize)
	var done bool
	//fastq input
	file := fileio.EasyOpen(input)
	defer file.Close()

	samRecords, _ := os.Create(output)
	defer samRecords.Close()

	header := NodesHeader(ref.Nodes)
	sam.WriteHeaderToFileHandle(samRecords, header)
	batchSize := 5
	batches := make([][]*fastq.Fastq, 0, batchSize)
	var batchID int = 0
	var wg sync.WaitGroup
	for fq, done = fastq.NextFastq(file); !done; fq, done = fastq.NextFastq(file) {
		groups = append(groups, fq)

		if len(groups) == groupSize {
			//send off group

			batches = append(batches, groups)
			groups = make([]*fastq.Fastq, 0, groupSize)

			if len(batches) == batchSize {
				log.Printf("Aligning reads...\n")
				alignBatchGroup(ref, batches, seedHash, seedLen, stepSize, 824, samRecords, batchID, &wg)
				batches = make([][]*fastq.Fastq, 0, batchSize)
			}
		}
	}
	//wg.Wait()
	//last case check length groups
	if len(batches) != 0 {
		if len(groups) != 0 {
			batches = append(batches, groups)
		}

		alignBatchGroup(ref, batches, seedHash, seedLen, stepSize, 824, samRecords, batchID, &wg)
	}
	//wg.Wait()
}*/

/*
func GenomeDiversitySimulator(dnaSequence string) {
	var tileSize int = 32
	var stepSize int = 31
	var numberOfReads int = 10000
	var readLength int = 150
	var mutations int = 0
	var workerWaiter, writerWaiter sync.WaitGroup
	var numWorkers int = 8

	log.Printf("Reading in the genome (simple graph)...\n")
	referenceGenome := Read(dnaSequence)
	log.Printf("Simulating %d reads...\n", numberOfReads)
	simReads := RandomReads(referenceGenome.Nodes, readLength, numberOfReads, mutations, true)
	//simReads := GenomeDiversity(genome.Nodes, readLength, numberOfReads)
	fastq.Write("genomeDiversity.fastq", simReads)
	log.Printf("Writing fastqs to file...\n")
	log.Printf("Reading in the genome (simple graph)...\n")
	genome := Read(dnaSequence)
	log.Printf("Indexing the genome...\n")
	tiles := IndexGenomeIntoMap(referenceGenome.Nodes, tileSize, stepSize)
	file, _ := os.Create("testdata/genomeDiversity.sam")
	defer file.Close()
	header := NodesHeader(genome.Nodes)
	sam.WriteHeaderToFileHandle(file, header)

	time.Sleep(10 * time.Second)

	log.Printf("Making fastq channel...\n")
	fastqPipe := make(chan *fastq.Fastq, 824)

	log.Printf("Making sam channel...\n")
	samPipe := make(chan *sam.SamAln, 824)

	log.Printf("Waiting for 10 seconds and then aligning reads...\n")
	time.Sleep(10 * time.Second)

	go fastq.ReadToChan("genomeDiversity.fastq", fastqPipe)

	log.Printf("Starting alignment worker...\n")
	workerWaiter.Add(numWorkers)
	start := time.Now()
	for i := 0; i < numWorkers; i++ {
		go gswWorker(genome, tiles, tileSize, stepSize, fastqPipe, samPipe, &workerWaiter)
	}

	writerWaiter.Add(1)

	go sam.SamChanToFile(samPipe, file, &writerWaiter)

	workerWaiter.Wait()
	close(samPipe)
	log.Printf("Aligners finished and channel closed\n")
	writerWaiter.Wait()
	log.Printf("Sam writer finished and we are all done\n")
	stop := time.Now()
	duration := stop.Sub(start)
	log.Printf("Aligned %d reads in %s (%.1f reads per second).\n", len(simReads), duration, float64(len(simReads))/duration.Seconds())
}*/
