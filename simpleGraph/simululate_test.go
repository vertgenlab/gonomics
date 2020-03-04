package simpleGraph

import (
	//"fmt"
	"github.com/vertgenlab/gonomics/cigar"
	"github.com/vertgenlab/gonomics/dna"
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

func TestSamToVcf(t *testing.T) {

	chrI := &Node{Name: "chrI", Seq: dna.StringToBases("ACAAAAAAAAAAAAA")}
	h := NodesHeader([]*Node{chrI})
	align1 := &sam.SamAln{QName: "", Flag: 0, RName: "chrI", Pos: 2, MapQ: 255, Cigar: cigar.FromString("1M5S"), RNext: "*", PNext: 0, TLen: 0, Seq: dna.StringToBases("GAAAAA"), Extra: "BZ:i:0"}
	align2 := &sam.SamAln{QName: "", Flag: 0, RName: "chrI", Pos: 2, MapQ: 255, Cigar: cigar.FromString("1M5S"), RNext: "*", PNext: 0, TLen: 0, Seq: dna.StringToBases("GAAAAA"), Extra: "BZ:i:0"}
	align3 := &sam.SamAln{QName: "", Flag: 0, RName: "chrI", Pos: 4, MapQ: 255, Cigar: cigar.FromString("4M"), RNext: "*", PNext: 0, TLen: 0, Seq: dna.StringToBases("AACT"), Extra: "BZ:i:0"}
	align4 := &sam.SamAln{QName: "", Flag: 0, RName: "chrI", Pos: 4, MapQ: 255, Cigar: cigar.FromString("4M"), RNext: "*", PNext: 0, TLen: 0, Seq: dna.StringToBases("AACT"), Extra: "BZ:i:0"}
	align5 := &sam.SamAln{QName: "", Flag: 0, RName: "chrI", Pos: 8, MapQ: 255, Cigar: cigar.FromString("1M7I"), RNext: "*", PNext: 0, TLen: 0, Seq: dna.StringToBases("AGATGAGT"), Extra: "BZ:i:0"}
	align6 := &sam.SamAln{QName: "", Flag: 0, RName: "chrI", Pos: 8, MapQ: 255, Cigar: cigar.FromString("1M7I"), RNext: "*", PNext: 0, TLen: 0, Seq: dna.StringToBases("AGATGAGT"), Extra: "BZ:i:0"}
	var testSam []*sam.SamAln = []*sam.SamAln{align1, align2, align3, align4, align5, align6}
	samTest := &sam.Sam{Header: h, Aln: testSam}
	sam.Write("testdata/samToVcfTest.sam", samTest)

	chrITest := &fasta.Fasta{Name: "chrI", Seq: dna.StringToBases("ACAAAAAAAAAAAAA")}
	votes := samToGenomeNotes("testdata/samToVcfTest.sam", chrITest, 0)
	log.Printf("Before: %s\n", dna.BasesToString(chrITest.Seq))
	fa, v := EditGenome("testdata/samToVcfTest.sam", chrITest, 0, votes)

	log.Printf("After:%s\n", dna.BasesToString(fa.Seq))
	vcf.PrintVcf(v)

}

func TestAlignPairedEnd(t *testing.T) {
	var tileSize int = 30
	var stepSize int = 24
	var numberOfReads int = 10000
	var readLength int = 150
	var mutations int = 0
	var workerWaiter, writerWaiter sync.WaitGroup
	var numWorkers int = 8

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

	header := NodesHeader(genome.Nodes)

	start := time.Now()
	go fastq.PairEndToChan("testdata/simReads_R1.fq", "testdata/simReads_R2.fq", fastqPipe)

	log.Printf("Starting alignment worker...\n")
	workerWaiter.Add(numWorkers)
	for i := 0; i < numWorkers; i++ {
		go gswPairEnd(genome, tiles, tileSize, stepSize, fastqPipe, samPipe, &workerWaiter)
	}
	writerWaiter.Add(1)
	go sam.SamChanPairToFile(samPipe, "/dev/stdout", header, &writerWaiter)
	workerWaiter.Wait()
	close(samPipe)
	log.Printf("Aligners finished and channel closed\n")
	writerWaiter.Wait()
	log.Printf("Sam writer finished and we are all done\n")
	stop := time.Now()
	duration := stop.Sub(start)
	os.Remove("samFile.sam")
	log.Printf("Aligned %d reads in %s (%.1f reads per second).\n", len(simReads)*2, duration, float64(len(simReads)*2)/duration.Seconds())
}

func TestVcfToGraph(t *testing.T) {
	var tileSize int = 32
	var stepSize int = 31
	var numberOfReads int = 20000
	var readLength int = 150
	//var mutations int = 0

	var workerWaiter, writerWaiter sync.WaitGroup
	var numWorkers int = 8
	log.Printf("Reading in the genome (fasta)...\n")
	ref := fasta.Read("testdata/gasAcu1.fa")
	fasta.AllToUpper(ref)
	//log.Printf("Simulating vcf...\n")
	//vcfFile := SimulateVcfGenomeWide(ref, mutations)
	//vcf.Sort(vcfFile)
	//vcf.PrintVcf(vcfFile)
	//log.Printf("Vcf to genome graph...\n")
	genome := Read("testdata/gasAcu1-4_snpIndels.gg")
	//VariantGraph(ref, vcfFile)
	//Write("testdata/vcfGraph.gg", genome)
	log.Printf("Simulating reads...\n")
	fa := Read("testdata/gasAcu1.fa")

	simReads := RandomReads(fa.Nodes, readLength, numberOfReads, 0)

	fastq.Write("testdata/simReads.fq", simReads)

	log.Printf("Indexing the genome...\n")
	tiles := IndexGenomeIntoMap(genome.Nodes, tileSize, stepSize)

	log.Printf("Making fastq channel...\n")
	fastqPipe := make(chan *fastq.Fastq, 824)

	log.Printf("Making sam channel...\n")
	samPipe := make(chan *sam.SamAln, 824)

	header := NodesHeader(genome.Nodes)

	start := time.Now()
	go fastq.ReadToChan("testdata/simReads.fq", fastqPipe)

	log.Printf("Starting alignment worker...\n")
	workerWaiter.Add(numWorkers)
	for i := 0; i < numWorkers; i++ {
		go gswWorker(genome, tiles, tileSize, stepSize, fastqPipe, samPipe, &workerWaiter)
	}
	writerWaiter.Add(1)
	go sam.SamChanToFile(samPipe, "/dev/stdout", header, &writerWaiter)
	workerWaiter.Wait()
	close(samPipe)
	log.Printf("Aligners finished and channel closed\n")
	writerWaiter.Wait()
	log.Printf("Sam writer finished and we are all done\n")
	stop := time.Now()
	duration := stop.Sub(start)
	//os.Remove("samFile.sam")
	log.Printf("Aligned %d reads in %s (%.1f reads per second).\n", len(simReads), duration, float64(len(simReads))/duration.Seconds())
	aln, _ := sam.Read("testdata/alignToGraph.sam")
	//for _, aligned := range aln.Aln{
	//	log.Printf("%s\n", LocalView(aligned, genome.Nodes))
	//}
	CheckAnswers(aln.Aln)

}
