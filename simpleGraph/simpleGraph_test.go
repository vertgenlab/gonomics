package simpleGraph

import (
	"flag"
	"fmt"
	//"github.com/vertgenlab/gonomics/common"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/fastq"
	"github.com/vertgenlab/gonomics/sam"
	"github.com/vertgenlab/gonomics/vcf"
	"log"
	//"math"
	"os"
	"runtime"
	"runtime/pprof"
	"sync"
	"testing"
	"time"
)

var seqOneA = dna.StringToBases("ACGTACGTCATCATCATTACTACTAC")
var seqOneB = dna.StringToBases("ACGTACGT")
var seqOneC = dna.StringToBases("ACGTACGTACGTT")
var readWriteTests = []struct {
	filename string // input
	//data *SimpleGraph
	data []*Node
}{
	//{"testdata/testOne.sg", []*Node{{0, "seqOneA", seqOneA, nil, nil}, {1, "seqOneB", seqOneB, nil, nil}, {2, "seqOneC", seqOneC, nil, nil}}},
}

/*
func TestPabBioGraph(t *testing.T) {
	var tileSize int = 32
	var stepSize int = 32
	var numberOfReads int = 10
	var readLength int = 150
	var mutations int = 0
	var workerWaiter, writerWaiter sync.WaitGroup
	var numWorkers int = 8
	var scoreMatrix = HumanChimpTwoScoreMatrix

	log.Printf("Reading in the genome (simple graph)...\n")
	genome, _ := Read("testdata/BRpbsv.gg")
	//genome, _ := Read("testdata/rabsBepaChrI.gg")

	log.Printf("Indexing the genome...\n")
	tiles := indexGenomeIntoMap(genome.Nodes, tileSize, stepSize)

	log.Printf("Making fastq channel...\n")
	fastqPipe := make(chan *fastq.FastqBig, 824)

	log.Printf("Making sam channel...\n")
	samPipe := make(chan *sam.SamAln, 824)

	log.Printf("Simulating reads...\n")
	simReads := RandomReads(genome, readLength, numberOfReads, mutations)
	fastq.Write("testdata/simReads.fq", simReads)

	header := NodesHeader(genome.Nodes)

	go fastq.ReadBigToChan("testdata/simReads.fq", fastqPipe)
	writerWaiter.Add(1)
	go sam.SamChanToFile(samPipe, "/dev/stdout", header, &writerWaiter)

	log.Printf("Starting alignment worker...\n")
	time.Sleep(5 * time.Second)

	f, err := os.Create("testdata/cpuprofile.data")
	common.ExitIfError(err)
	defer f.Close()
	err = pprof.StartCPUProfile(f)
	common.ExitIfError(err)

	start := time.Now()
	workerWaiter.Add(numWorkers)
	for i := 0; i < numWorkers; i++ {
		go gswWorkerMemPool(genome, tiles, tileSize, stepSize, scoreMatrix, fastqPipe, samPipe, &workerWaiter)
	}

	workerWaiter.Wait()
	stop := time.Now()
	pprof.StopCPUProfile()
	close(samPipe)

	log.Printf("Aligners finished and channel closed\n")
	writerWaiter.Wait()
	log.Printf("Sam writer finished and we are all done\n")
	duration := stop.Sub(start)
	log.Printf("Aligned %d reads in %s (%.1f reads per second).\n", len(simReads), duration, float64(len(simReads))/duration.Seconds())
	//strictCheckFile("testPabBioGraph.sam", genome)
	//log.Printf("Passed alignment check!!!\n")
}*/

func TestWorkerWithWriting(t *testing.T) {
	var tileSize int = 32
	var stepSize int = 31
	var numberOfReads int = 10
	var readLength int = 150
	var mutations int = 0
	var workerWaiter, writerWaiter sync.WaitGroup
	var numWorkers int = 3
	genome, _ := Read("testdata/bigGenome.sg")
	log.Printf("Reading in the genome (simple graph)...\n")
	//genome := Read("testdata/rabs_chr8.gg")
	//fa, _ := Read("testdata/chrI.fa")
	log.Printf("Indexing the genome...\n")
	log.Printf("Making fastq channel...\n")
	fastqPipe := make(chan *fastq.Fastq, 824)

	log.Printf("Making sam channel...\n")
	samPipe := make(chan *sam.SamAln, 824)

	log.Printf("Simulating reads...\n")
	simReads := RandomReads(genome, readLength, numberOfReads, mutations)
	os.Remove("testdata/simReads.fq")
	fastq.Write("testdata/simReads.fq", simReads)
	//genome, chrSize := Read("testdata/rabsBepaChrI.gg")

	//header := sam.ChromInfoMapSamHeader(chrSize)

	//log.Printf("%v\n", header)
	tiles := IndexGenomeIntoMap(genome.Nodes, tileSize, stepSize)
	log.Printf("Finished Indexing Genome...\n")
	start := time.Now()
	go fastq.ReadToChan("testdata/simReads.fq", fastqPipe)

	log.Printf("Starting alignment worker...\n")
	workerWaiter.Add(numWorkers)
	for i := 0; i < numWorkers; i++ {
		go gswWorker(genome, tiles, tileSize, stepSize, fastqPipe, samPipe, &workerWaiter)
	}
	writerWaiter.Add(1)
	go SamChanView(samPipe, genome, &writerWaiter)
	///dev/stdout
	//go sam.SamChanToFile(samPipe, "/dev/stdout", header, &writerWaiter)
	workerWaiter.Wait()
	close(samPipe)
	log.Printf("Aligners finished and channel closed\n")
	writerWaiter.Wait()
	log.Printf("Sam writer finished and we are all done\n")
	stop := time.Now()
	duration := stop.Sub(start)
	//os.Remove("testdata/sim.sam")
	log.Printf("Aligned %d reads in %s (%.1f reads per second).\n", len(simReads), duration, float64(len(simReads))/duration.Seconds())
}

/*func TestWorkerWithTiming(t *testing.T) {
	var tileSize int = 32
	var stepSize int = tileSize - 1
	var numberOfReads int = 100
	var readLength int = 150
	var mutations int = 2
	var dummyWaiter sync.WaitGroup

	log.Printf("Reading in the genome (simple graph)...\n")
	fa, _ := Read("testdata/bigGenome.sg")

	log.Printf("Making fastq channel...\n")
	fastqPipe := make(chan *fastq.Fastq, 824)

	log.Printf("Making sam channel...\n")
	samPipe := make(chan *sam.SamAln, 824)

	log.Printf("Simulating reads...\n")
	simReads := RandomReads(fa.Nodes, readLength, numberOfReads, mutations)
	//dummyWaiter.Add(len(simReads))
	alignments := make([]*sam.SamAln, numberOfReads)
	//TODO: this file does not exist in testdata
	genome, _ := Read("testdata/rabsToGasAcu1.gg")
	log.Printf("Indexing the genome...\n")

	tiles := IndexGenomeIntoMap(genome.Nodes, tileSize, stepSize)

	log.Printf("Starting alignment worker...\n")
	dummyWaiter.Add(1)
	go gswWorker(genome, tiles, tileSize, stepSize, fastqPipe, samPipe, &dummyWaiter)

	log.Printf("Waiting for 5 seconds and then aligning reads...\n")
	time.Sleep(5 * time.Second)

	start := time.Now()
	for i := 0; i < numberOfReads; i++ {
		fastqPipe <- simReads[i]
	}
	for j := 0; j < numberOfReads; j++ {
		alignments[j] = <-samPipe
	}
	dummyWaiter.Wait()
	stop := time.Now()
	duration := stop.Sub(start)
	log.Printf("Aligned %d reads in %s (%.1f reads per second).\n", len(alignments), duration, float64(len(alignments))/duration.Seconds())
	//for k := 0; k < len(alignments); k++ {
	//fmt.Printf("%s\n", LocalView(alignments[k], genome.Nodes))
	//	fmt.Printf("%s\n", ViewGraphAignment(alignments[k], genome))
	//}
	//CheckAnswers(alignments, genome)
}*/

func TestHippoAln(t *testing.T) {
	var hippo *fastq.Fastq = &fastq.Fastq{Name: "hippoOne", Seq: dna.StringToBases("TGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGAGTGATTTGAAGGTACATGGAATACCACCACGGGAGCAAAGC"), Qual: []rune("JJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJ")}
	var tileSize int = 32
	var stepSize int = tileSize - 1
	var alignment *sam.SamAln = nil
	var dummyWaiter sync.WaitGroup

	log.Printf("Reading in the genome (simple graph)...\n")
	genome, _ := Read("testdata/bigGenome.sg")

	log.Printf("Indexing the genome...\n")
	tiles := IndexGenomeIntoMap(genome.Nodes, tileSize, stepSize)

	log.Printf("Making fastq channel...\n")
	fastqPipe := make(chan *fastq.Fastq, 1)

	log.Printf("Making sam channel...\n")
	samPipe := make(chan *sam.SamAln, 1)

	log.Printf("Starting alignment worker...\n")
	go gswWorker(genome, tiles, tileSize, stepSize, fastqPipe, samPipe, &dummyWaiter)

	log.Printf("Waiting for 5 seconds and then aligning read...\n")
	time.Sleep(5 * time.Second)

	start := time.Now()
	fastqPipe <- hippo
	alignment = <-samPipe
	end := time.Now()
	duration := end.Sub(start)
	log.Printf("duration:%s\t%s\n", duration, dna.BasesToString(alignment.Seq))
}

func TestAligning(t *testing.T) {
	var tileSize int = 32
	var stepSize int = tileSize - 1
	var readLength int = 150
	var numberOfReads int = 5
	var mutations int = 0
	log.Printf("Reading in the genome (simple graph)...\n")
	genome, _ := Read("testdata/bigGenome.sg")

	log.Printf("Indexing the genome...\n")
	tiles := IndexGenomeIntoMap(genome.Nodes, tileSize, stepSize)

	log.Printf("Simulating reads...\n")
	simReads := RandomReads(genome, readLength, numberOfReads, mutations)

	log.Printf("Aligning reads...\n")
	var cpuprofile = flag.String("cpuprofile", "", "write cpu profile to `file`")
	var memprofile = flag.String("memprofile", "", "write memory profile to `file`")

	flag.Parse()
	if *cpuprofile != "" {
		f, err := os.Create(*cpuprofile)
		if err != nil {
			log.Fatal("could not create CPU profile: ", err)
		}
		defer f.Close()
		if err := pprof.StartCPUProfile(f); err != nil {
			log.Fatal("could not start CPU profile: ", err)
		}
		defer pprof.StopCPUProfile()
	}

	start := time.Now()
	for i := 0; i < len(simReads); i++ {
		before := time.Now()
		wrapNoChanMap(genome, simReads[i], tiles, tileSize, stepSize)
		after := time.Now()
		timeForOne := after.Sub(before)
		log.Printf("It took %s to map the last read\n", timeForOne)
	}
	stop := time.Now()
	elapsed := stop.Sub(start)
	//code block for program profiling
	if *memprofile != "" {
		f, err := os.Create(*memprofile)
		if err != nil {
			log.Fatal("could not create memory profile: ", err)
		}
		defer f.Close()
		runtime.GC() // get up-to-date statistics
		if err := pprof.WriteHeapProfile(f); err != nil {
			log.Fatal("could not write memory profile: ", err)
		}
	}

	log.Printf("It took %s to map %d reads\n", elapsed, numberOfReads)
}

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

func TestVcfGraph(t *testing.T) {
	//smallFasta := fasta.Fasta{Name: "chr1", Seq: dna.StringToBases("ATCGA")}
	smallFasta := fasta.Fasta{Name: "chr1", Seq: dna.StringToBases("ATCGA")}
	//smallFasta := fasta.Fasta{Name: "chr1", Seq: dna.StringToBases("ATCGA")}
	fmt.Printf("Reference sequence is: %s\n", dna.BasesToString(smallFasta.Seq))
	var vcfTest []*vcf.Vcf
	//vcfTest = append(vcfTest, &vcf.Vcf{Chr: "chr1", Pos: 3, Id: ".", Ref: "T", Alt: "TA", Qual: 0, Filter: "PASS", Info: "", Format: "SVTYPE=INS"})
	vcfTest = append(vcfTest, &vcf.Vcf{Chr: "chr1", Pos: 3, Id: ".", Ref: "C", Alt: "T", Qual: 0, Filter: "PASS", Info: "", Format: "SVTYPE=SNP"})
	//vcfTest = append(vcfTest, &vcf.Vcf{Chr: "chr1", Pos: 3, Id: ".", Ref: "C", Alt: "CTAG", Qual: 0, Filter: "PASS", Info: "", Format: "SVTYPE=INS"})
	//vcfTest = append(vcfTest, &vcf.Vcf{Chr: "chr1", Pos: 2, Id: ".", Ref: "TCG", Alt: "T", Qual: 0, Filter: "PASS", Info: "", Format: "SVTYPE=DEL"})
	//vcfTest = append(vcfTest, &vcf.Vcf{Chr: "chr1", Pos: 14, Id: ".", Ref: "G", Alt: "GAAAAAAAAAA", Qual: 0, Filter: "PASS", Info: "", Format: "SVTYPE=INS"})
	//vcfTest = append(vcfTest, &vcf.Vcf{Chr: "chr1", Pos: 19, Id: ".", Ref: "A", Alt: "ATTTTT", Qual: 0, Filter: "PASS", Info: "", Format: "SVTYPE=INS"})
	//var sg *SimpleGraph = NewGraph()

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
	newSg, _ := Read("testdata/vcfGraphTest.gg")

	PrintGraph(newSg)
	vcf.Write("anotherTesting.vcf", vcfTest)
}

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
}

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
func TestGraphTraversal(t *testing.T) {
	gg := NewGraph()

	var seqOne = dna.StringToBases("ATG")
	var seqTwo = dna.StringToBases("GGA")
	var seqThree = dna.StringToBases("CCC")
	var seqFour = dna.StringToBases("TTG")
	var seqFive = dna.StringToBases("TGT")
	m, trace := swMatrixSetup(10000)
	testFastq := fastq.Fastq{Name: "TestSeq", Seq: dna.StringToBases("GGA"), Qual: []rune("JJJ")}
	bestScore, cigs, refStart, refEnd, queryStart, queryEnd := LeftLocal(dna.StringToBases("AAATG"), dna.StringToBases("TG"), HumanChimpTwoScoreMatrix, -600, m, trace)
	fmt.Printf("bestScore: %d, cigs: %v, refStart: %d, refEnd: %d, queryStart: %d, queryEnd: %d", bestScore, cigs, refStart, refEnd, queryStart, queryEnd)
	nA := Node{0, "A", seqOne, nil, nil}
	nB := Node{1, "B", seqTwo, nil, nil}
	nC := Node{3, "C", seqThree, nil, nil}
	nD := Node{4, "D", seqFour, nil, nil}
	nE := Node{5, "E", seqFive, nil, nil}

	AddNode(gg, &nA)
	AddNode(gg, &nB)
	AddNode(gg, &nC)
	AddNode(gg, &nD)
	AddNode(gg, &nE)

	//AddEdge(gg, &nA, &nB, 1)
	//AddEdge(gg, &nB, &nC, 1)

	AddEdge(&nA, &nB, 1)
	AddEdge(&nA, &nC, 1)
	AddEdge(&nA, &nE, 1)
	AddEdge(&nB, &nD, 1)
	AddEdge(&nC, &nD, 1)
	AddEdge(&nE, &nD, 1)
	var query []dna.Base
	var path []uint32
	//fmt.Printf("Start node: %s\n", dna.BasesToString(gg.Nodes[3].Seq))
	path, query = ReverseGraphTraversal(gg.Nodes[1], query, path, 1, 9)
	reversePath(path)
	log.Printf("Reverse print path: %s", PathToString(path, gg))
	GraphTraversalFwd(gg, gg.Nodes[1], query, path, 1, 9)

	var mappedRead *sam.SamAln = &sam.SamAln{QName: testFastq.Name, Flag: 0, RName: "", Pos: 0, MapQ: 255, RNext: "*", PNext: 0, TLen: 0, Seq: []dna.Base{}, Qual: "", Extra: ""}
	var alignment []*cigar.Cigar
	var score int64 = 0
	var tStart, qStart int
	var testCig []*cigar.Cigar
	//var cig []*cigar.Cigar
	log.Println("Starting the forward alignment...")
	AlignTraversalFwd(gg.Nodes[1], query, 0, path, 0, testFastq.Seq, m, trace)
	fmt.Printf("traveral final: %v %d\n", testCig, score)
	score = 0
	alignment, score, tStart, qStart, path = AlignReverseGraphTraversal(gg.Nodes[1], query, 0, path, 6, testFastq.Seq, m, trace)

	fmt.Printf("Reverse score: %d query start: %d\n", score, qStart)
	mappedRead.Pos = int64(tStart)
	mappedRead.Cigar = alignment
	log.Printf("%s\n", sam.SamAlnToString(mappedRead))
	PrintGraph(gg)

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
