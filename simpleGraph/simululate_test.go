package simpleGraph

import (
	//"fmt"
	//"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/fastq"
	"github.com/vertgenlab/gonomics/sam"
	"github.com/vertgenlab/gonomics/vcf"
	"log"
	//"math"
	"os"
	"sync"
	"testing"
	"time"
)

func TestAlignPairedEnd(t *testing.T) {
	var tileSize int = 32
	var stepSize int = 31
	var numberOfReads int = 10000
	var readLength int = 150
	var mutations int = 0
	var workerWaiter, writerWaiter sync.WaitGroup
	var numWorkers int = 4

	log.Printf("Reading in the genome (Simple Graph)...\n")
	genome := Read("testdata/gasAcu1.fa")
	log.Printf("Simulating reads...\n")

	simReads := PairedEndRandomReads(genome.Nodes, readLength, numberOfReads, mutations)
	log.Printf("length of simulated paired end reads: %d\n", len(simReads))
	fastq.WritePair("testdata/simReads_R1.fq", "testdata/simReads_R2.fq", simReads)
	log.Printf("Making fastq channel...\n")
	fastqPipe := make(chan *fastq.PairedEnd, 824)

	log.Printf("Indexing the genome...\n")
	tiles := IndexGenomeIntoMap(genome.Nodes, tileSize, stepSize)

	log.Printf("Making sam channel...\n")
	samPipe := make(chan *sam.PairedSamAln, 824)

	file, _ := os.Create("/dev/stdout")
	defer file.Close()
	header := NodesHeader(genome.Nodes)
	sam.WriteHeaderToFileHandle(file, header)

	start := time.Now()
	go fastq.PairEndToChan("testdata/simReads_R1.fq", "testdata/simReads_R2.fq", fastqPipe)

	log.Printf("Starting alignment worker...\n")
	workerWaiter.Add(numWorkers)
	for i := 0; i < numWorkers; i++ {
		go gswPairEnd(genome, tiles, tileSize, stepSize, fastqPipe, samPipe, &workerWaiter)
	}
	writerWaiter.Add(1)
	go sam.SamChanPairToFile(samPipe, file, &writerWaiter)
	workerWaiter.Wait()
	close(samPipe)
	log.Printf("Aligners finished and channel closed\n")
	writerWaiter.Wait()
	log.Printf("Sam writer finished and we are all done\n")
	stop := time.Now()
	duration := stop.Sub(start)
	os.Remove("samFile.sam")
	log.Printf("Aligned %d reads in %s (%.1f reads per second).\n", len(simReads), duration, float64(len(simReads))/duration.Seconds())
}

func TestVcfToGraph(t *testing.T) {
	var tileSize int = 32
	var stepSize int = 31
	var numberOfReads int = 10000
	var readLength int = 150
	var mutations int = 1000
	var workerWaiter, writerWaiter sync.WaitGroup
	var numWorkers int = 4
	log.Printf("Reading in the genome (fasta)...\n")
	ref := fasta.Read("testdata/gasAcu1.fa")
	log.Printf("Simulating vcf...\n")
	vcfFile := SimulateVcfGenomeWide(ref, mutations)
	vcf.PrintVcf(vcfFile)
	log.Printf("Vcf to genome graph...\n")
	genome := FaToGenomeGraph(ref, vcfFile)
	//Write("testdata/vcfGraph.gg", genome)
	log.Printf("Simulating reads...\n")
	fa := Read("testdata/gasAcu1.fa")
	simReads := RandomReads(fa.Nodes, readLength, numberOfReads, mutations)

	fastq.Write("testdata/simReads.fq", simReads)

	log.Printf("Indexing the genome...\n")
	tiles := IndexGenomeIntoMap(genome.Nodes, tileSize, stepSize)

	log.Printf("Making fastq channel...\n")
	fastqPipe := make(chan *fastq.Fastq, 824)

	log.Printf("Making sam channel...\n")
	samPipe := make(chan *sam.SamAln, 824)

	file, _ := os.Create("samFile.sam")
	defer file.Close()
	header := NodesHeader(genome.Nodes)
	sam.WriteHeaderToFileHandle(file, header)

	start := time.Now()
	go fastq.ReadToChan("testdata/simReads.fq", fastqPipe)

	log.Printf("Starting alignment worker...\n")
	workerWaiter.Add(numWorkers)
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
	os.Remove("samFile.sam")
	log.Printf("Aligned %d reads in %s (%.1f reads per second).\n", len(simReads), duration, float64(len(simReads))/duration.Seconds())
}
