package simpleGraph

import (
	"github.com/vertgenlab/gonomics/fastq"
	"github.com/vertgenlab/gonomics/giraf"
	"log"
	"os"
	"sync"
	"testing"
	"time"
)

func TestWorkerWithWriting(t *testing.T) {
	var output string = "testdata/pairedTest.giraf"
	var tileSize int = 32
	var stepSize int = 32
	var numberOfReads int = 100
	var readLength int = 150
	var mutations int = 0
	var workerWaiter, writerWaiter sync.WaitGroup
	var numWorkers int = 8
	var scoreMatrix = HumanChimpTwoScoreMatrix
	genome := Read("testdata/bigGenome.sg")
	log.Printf("Reading in the genome (simple graph)...\n")
	log.Printf("Indexing the genome...\n")
	log.Printf("Making fastq channel...\n")
	fastqPipe := make(chan *fastq.PairedEndBig, 824)

	log.Printf("Making sam channel...\n")
	samPipe := make(chan *giraf.GirafPair, 824)

	log.Printf("Simulating reads...\n")
	simReads := RandomPairedReads(genome, readLength, numberOfReads, mutations)
	os.Remove("testdata/simReads_R1.fq")
	os.Remove("testdata/simReads_R2.fq")
	fastq.WritePair("testdata/simReads_R1.fq", "testdata/simReads_R2.fq", simReads)
	tiles := IndexGenomeIntoMap(genome.Nodes, tileSize, stepSize)
	go fastq.ReadPairBigToChan("testdata/simReads_R1.fq", "testdata/simReads_R2.fq", fastqPipe)
	log.Printf("Finished Indexing Genome...\n")
	start := time.Now()

	log.Printf("Starting alignment worker...\n")
	workerWaiter.Add(numWorkers)
	for i := 0; i < numWorkers; i++ {
		go RoutineFqPairToGiraf(genome, tiles, tileSize, stepSize, scoreMatrix, fastqPipe, samPipe, &workerWaiter)
	}
	go giraf.GirafPairChanToFile(output, samPipe, &writerWaiter)
	writerWaiter.Add(1)
	workerWaiter.Wait()
	close(samPipe)
	log.Printf("Aligners finished and channel closed\n")
	writerWaiter.Wait()
	log.Printf("Sam writer finished and we are all done\n")
	stop := time.Now()
	duration := stop.Sub(start)
	log.Printf("Aligned %d reads in %s (%.1f reads per second).\n", len(simReads)*2, duration, float64(len(simReads)*2)/duration.Seconds())

}

/*
func TestHippoAln(t *testing.T) {
	var hippo *fastq.Fastq = &fastq.Fastq{Name: "hippoOne", Seq: dna.StringToBases("TGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGAGTGATTTGAAGGTACATGGAATACCACCACGGGAGCAAAGC"), Qual: fastq.ToQualUint8([]rune("JJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJ"))}
	var tileSize int = 32
	var stepSize int = tileSize - 1
	var alignment *sam.SamAln = nil
	var dummyWaiter sync.WaitGroup
	var scoreMatrix = HumanChimpTwoScoreMatrix
	log.Printf("Reading in the genome (simple graph)...\n")
	genome := Read("testdata/bigGenome.sg")

	log.Printf("Indexing the genome...\n")
	tiles := IndexGenomeIntoMap(genome.Nodes, tileSize, stepSize)

	log.Printf("Making fastq channel...\n")
	fastqPipe := make(chan *fastq.FastqBig, 1)

	log.Printf("Making sam channel...\n")
	samPipe := make(chan *giraf.Giraf, 1)

	log.Printf("Starting alignment worker...\n")
	go RoutineFqToGiraf(genome, tiles, scoreMatrix, tileSize, stepSize, fastqPipe, samPipe, &dummyWaiter)

	log.Printf("Waiting for 5 seconds and then aligning read...\n")
	time.Sleep(5 * time.Second)

	start := time.Now()
	fastqPipe <- hippo
	alignment = <-samPipe
	end := time.Now()
	duration := end.Sub(start)
	log.Printf("duration:%s\t%s\n", duration, dna.BasesToString(alignment.Seq))
}*/

/*func TestReadsWithTiming(t *testing.T) {
	//var hippo *fastq.Fastq = &fastq.Fastq{Name: "hippoOne", Seq: dna.StringToBases("ACCTTTTTCTTGTTGTATTTAAAGACAAATGATTTGATTTTATATAGCCAAATGGTTTTCAACGCTAGCAGTGTTTGGTGGCAACTCAGTTTCACCCACGTCTGTTCCAACTAACATGCAATATGTTTCCTGTAATCTGCAGCACGCTTT"), Qual: []rune("JJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJ")}
	var tileSize int = 32
	var stepSize int = 31
	var numberOfReads int = 1
	var readLength int = 150
	var mutations int = 2
	var dummyWaiter sync.WaitGroup

	var fastestRead, slowestRead *fastq.Fastq = nil, nil
	var fastestTime, slowestTime float64 = math.MaxFloat64, 0

	log.Printf("Reading in the genome (simple graph)...\n")
	fa, _ := Read("testdata/bigGenome.sg")
	//TODO: this file does not exist in testdata
	genome, _ := Read("testdata/rabsBepa.gg")
	log.Printf("Simulating reads...\n")
	simReads := RandomReads(fa.Nodes, readLength, numberOfReads, mutations)

	log.Printf("Indexing the genome...\n")
	tiles := IndexGenomeIntoMap(genome.Nodes, tileSize, stepSize)

	log.Printf("Making fastq channel...\n")
	fastqPipe := make(chan *fastq.Fastq, numberOfReads)

	log.Printf("Making sam channel...\n")
	samPipe := make(chan *sam.SamAln, numberOfReads)

	alignments := make([]*sam.SamAln, numberOfReads)

	log.Printf("Starting alignment worker...\n")

	go gswWorker(genome, tiles, tileSize, stepSize, fastqPipe, samPipe, &dummyWaiter)

	log.Printf("Waiting for 5 seconds and then aligning reads...\n")
	time.Sleep(5 * time.Second)

	for i := 0; i < numberOfReads; i++ {

		start := time.Now()
		fastqPipe <- simReads[i]
		alignments[i] = <-samPipe
		stop := time.Now()
		duration := stop.Sub(start).Seconds()
		if duration > slowestTime {
			slowestTime = duration
			slowestRead = simReads[i]
		} else if duration < fastestTime {
			fastestTime = duration
			fastestRead = simReads[i]
		}
	}
	log.Printf("Fastest read was (%.4f):\n%s\nSlowest reads was (%.4f):\n%s\n", fastestTime, dna.BasesToString(fastestRead.Seq), slowestTime, dna.BasesToString(slowestRead.Seq))

	//CheckAnswers(alignments, genome)
}*/

// This looks like a good test, but it is failing when trying to
// create the graph from the smallFasta and vcfTest
/*func TestVcfGraph(t *testing.T) {
	//smallFasta := fasta.Fasta{Name: "chr1", Seq: dna.StringToBases("ATCGA")}
	smallFasta := fasta.Fasta{Name: "chr1", Seq: dna.StringToBases("ATCGA")}
	//smallFasta := fasta.Fasta{Name: "chr1", Seq: dna.StringToBases("ATCGA")}
	fmt.Printf("Reference sequence is: %s\n", dna.BasesToString(smallFasta.Seq))
	var vcfTest []*vcf.Vcf
	//vcfTest = append(vcfTest, &vcf.Vcf{Chr: "chr1", Pos: 3, Id: ".", Ref: "T", Alt: "TA", Qual: 0, Filter: "PASS", Info: "", Format: "SVTYPE=INS"})
	vcfTest = append(vcfTest, &vcf.Vcf{Chr: "chr1", Pos: 3, Id: ".", Ref: "C", Alt: "T", Qual: 0, Filter: "PASS", Info: "", Format: "SVTYPE=SNP"})
	sg := vChrGraph(NewGraph(), &smallFasta, vcfTest)

	for n := 0; n < len(sg.Nodes); n++ {
		fmt.Printf("Printing nodes: %d, seq=%s, numOfEdges=%d ", sg.Nodes[n].Id, dna.BasesToString(sg.Nodes[n].Seq), len(sg.Nodes[n].Next))
		for e := 0; e < len(sg.Nodes[n].Next); e++ {
			fmt.Printf("-> %v, weight=%v", dna.BasesToString(sg.Nodes[n].Next[e].Dest.Seq), sg.Nodes[n].Next[e].Prob)
		}
		fmt.Println("")
	}
	PrintGraph(sg)
	Write("testdata/vcfGraphTest.gg", sg)
	newSg := Read("testdata/vcfGraphTest.gg")

	PrintGraph(newSg)
	vcf.Write("anotherTesting.vcf", vcfTest)
}*/

/*
func BenchmarkGoRoutinesMap(b *testing.B) {
	var tileSize int = 32
	var stepSize int = tileSize - 1
	var readLength int = 150
	var numberOfReads int = 5
	var mutations int = 0

	genome, _ := Read("testdata/bigGenome.sg")
	tiles := IndexGenomeIntoMap(genome.Nodes, tileSize, stepSize)
	simReads := RandomReads(genome, readLength, numberOfReads, mutations)
	b.ResetTimer()

	for n := 0; n < b.N; n++ {
		for i := 0; i < len(simReads); i++ {
			wrapNoChanMap(genome, simReads[i], tiles, tileSize, stepSize)
		}
	}
}*/

//TODO: slices not working right now
/*func BenchmarkGoRoutinesSlice(b *testing.B) {
	var tileSize int = 12
	var stepSize int = tileSize - 1
	var readLength int = 150
	var numberOfReads int = 50
	var mutations int = 0

	genome, _ := Read("testdata/bigGenome.sg")
	tiles := IndexGenomeIntoSlice(genome.Nodes, tileSize, stepSize)
	simReads := RandomReads(genome.Nodes, readLength, numberOfReads, mutations)
	b.ResetTimer()

	for n := 0; n < b.N; n++ {
		for i := 0; i < len(simReads); i++ {
			wrapNoChan(genome, simReads[i], tiles, tileSize)
		}
	}
}*/

/*
func TestFaFormat(t *testing.T) {

	var tileSize int = 32
	var numberOfReads int = 100
	var readLength int = 150
	var mutations int = 0
	var numWorkers int = 8
	log.Printf("Reading in the genome (simple graph)...\n")
	genome, chrSize := Read("testdata/gasAcu1.fa")
	simReads := RandomPairedReads(genome.Nodes, readLength, numberOfReads, mutations)
	fastq.WritePair("testdata/simReads_R1.fq", "testdata/simReads_R2.fq", simReads)
	header := sam.ChromInfoMapSamHeader(chrSize)
	GswAlignFaFormat(genome, "testdata/simReads_R1.fq", "testdata/simReads_R2.fq", "testdata/format.sam", numWorkers, tileSize, header)
}*/

/*
func TestNewSeed(t *testing.T) {
	var tileSize int = 32
	var stepSize int = tileSize - 1
	var numberOfReads int = 30000
	var readLength int = 150
	var mutations int = 0
	var workerWaiter, writerWaiter sync.WaitGroup
	var numWorkers int = 4

	log.Printf("Reading in the genome (simple graph)...\n")
	genome := Read("testdata/gasAcu1.fa")

	log.Printf("Indexing the genome...\n")
	tiles := IndexGenomeIntoMap(genome.Nodes, tileSize, stepSize)

	log.Printf("Making fastq channel...\n")
	fastqPipe := make(chan *fastq.Fastq, 824)

	log.Printf("Making sam channel...\n")
	samPipe := make(chan *sam.SamAln, 824)

	log.Printf("Simulating reads...\n")
	simReads := RandomReads(genome.Nodes, readLength, numberOfReads, mutations)
	fastq.Write("testdata/simReads.fq", simReads)
	start := time.Now()
	go fastq.ReadToChan("testdata/simReads.fq", fastqPipe)

	log.Printf("Starting alignment worker...\n")
	workerWaiter.Add(numWorkers)
	for i := 0; i < numWorkers; i++ {
		go SeedTestGswWorker(genome, tiles, tileSize, stepSize, fastqPipe, samPipe, &workerWaiter)
	}

	writerWaiter.Add(1)
	go sam.SamChanToFile(samPipe, "/dev/stdout", &writerWaiter)

	workerWaiter.Wait()
	close(samPipe)
	log.Printf("Aligners finished and channel closed\n")
	writerWaiter.Wait()
	log.Printf("Sam writer finished and we are all done\n")
	stop := time.Now()
	duration := stop.Sub(start)
	log.Printf("Aligned %d reads in %s (%.1f reads per second).\n", len(simReads), duration, float64(len(simReads))/duration.Seconds())
}*/

/*
func BenchmarkAligningNoGoroutines(b *testing.B) {
	//var mappedReads []*sam.SamAln = make([]*sam.SamAln, numberOfReads)
	numberOfReads = 10
	genome := Read("testdata/bigGenome.sg")
	//tiles := indexGenome(genome.Nodes, tileSize)
	tiles := IndexGenomeDev(genome.Nodes, tileSize, stepSize)
	simReads := RandomReads(genome.Nodes, readLength, numberOfReads, mutations)
	fastq.Write("testdata/fakeReads.fastq", simReads)
	//var seeds []Seed = make([]Seed, 256)
	m, trace := swMatrixSetup(10000)
	//var seeds []Seed = make([]Seed, 256)
	b.ResetTimer()

	//c := make(chan *sam.SamAln)
	var mappedRead *sam.SamAln
	for n := 0; n < b.N; n++ {
		for i := 0; i < len(simReads);i++ {
			mappedRead = GraphSmithWaterman(genome, simReads[i], tiles, tileSize, m, trace)
			log.Printf("%s\n", sam.SamAlnToString(mappedRead))
		}

		//devGSWsBatch(genome, "testdata/fakeReads.fastq", tiles, tileSize, m, trace, 400)
		//noLimitGSW(genome, "testdata/fakeReads.fastq", tiles, tileSize, m, trace)
	}
}*/

/*
func TestRead(t *testing.T) {
	for _, test := range readWriteTests {
		actual := Read(test.filename)

		if !AllAreEqual(test.data.Nodes, actual.Nodes) {
			t.Errorf("The %s file was not read correctly.", test.filename)
		}
	}
}

func TestWriteAndRead(t *testing.T) {
	//var actual []*Node
	for _, test := range readWriteTests {
		tempFile := test.filename + ".tmp"
		Write(tempFile, test.data)
		actual := Read(tempFile)
		if !AllAreEqual(test.data, actual.Nodes) {
			t.Errorf("The %s file was not read correctly.", test.filename)
		}
		err := os.Remove(tempFile)
		if err != nil {
			t.Errorf("Deleting temp file %s gave an error.", tempFile)
		}
	}
}*/
