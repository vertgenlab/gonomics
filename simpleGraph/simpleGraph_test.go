package simpleGraph

import (
	"github.com/vertgenlab/gonomics/fastq"
	"github.com/vertgenlab/gonomics/giraf"

	"log"
	"os"
	"sync"
	"testing"
	//"net/http"
	//"net/http/pprof"
	"flag"
	"runtime"
	"runtime/pprof"
	"time"
)

var cpuprofile = flag.String("cpuprofile", "", "write cpu profile to `file`")
var memprofile = flag.String("memprofile", "", "write memory profile to `file`")
/*
//test code in master right now
func TestWorkerWithWritingMASTER(t *testing.T) {
	var output string = "testdata/pairedTest.giraf"
	var tileSize int = 32
	var stepSize int = 32
	var numberOfReads int = 5000
	var readLength int = 150
	var mutations int = 2
	var workerWaiter, writerWaiter sync.WaitGroup
	var numWorkers int = 8
	var scoreMatrix = HumanChimpTwoScoreMatrix
	if *cpuprofile != "" {
		f, err := os.Create(*cpuprofile)
		if err != nil {
			log.Fatal("could not create CPU profile: ", err)
		}
		defer f.Close() // error handling omitted for example
		if err := pprof.StartCPUProfile(f); err != nil {
			log.Fatal("could not start CPU profile: ", err)
		}
		defer pprof.StopCPUProfile()
	}
	genome := Read("testdata/rabsTHREEspine.fa")
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
	if *memprofile != "" {
		f, err := os.Create(*memprofile)
		if err != nil {
			log.Fatal("could not create memory profile: ", err)
		}
		defer f.Close() // error handling omitted for example
		runtime.GC()    // get up-to-date statistics
		if err := pprof.WriteHeapProfile(f); err != nil {
			log.Fatal("could not write memory profile: ", err)
		}
	}
}*/
/*()
func TestNewGraph(t *testing.T) {

	var output string = "testdata/test.giraf"
	var tileSize int = 32
	var stepSize int = 32
	var numberOfReads int = 50000
	var readLength int = 150
	var mutations int = 0
	var workerWaiter, writerWaiter sync.WaitGroup
	var numWorkers int = 8
	var scoreMatrix = HumanChimpTwoScoreMatrix
	genome := SimplyRead("testdata/rabsTHREEspine.fa")

	log.Printf("Simulating reads...\n")
	simReads := RandomPairedReads(&genome, readLength, numberOfReads, mutations)
	os.Remove("testdata/simReads_R1.fq")
	os.Remove("testdata/simReads_R2.fq")
	fastq.WritePair("testdata/simReads_R1.fq", "testdata/simReads_R2.fq", simReads)

	log.Printf("Reading in the genome (simple graph)...\n")
	log.Printf("Indexing the genome...\n")
	log.Printf("Making fastq channel...\n")
	fastqPipe := make(chan *fastq.PairedEndBig, 2408)

	log.Printf("Making giraf channel...\n")
	girafPipe := make(chan *giraf.GirafPair, 2408)

	tiles := IndexGenomeIntoMap(genome.Nodes, tileSize, stepSize)
	go readFastqGsw("testdata/simReads_R1.fq", "testdata/simReads_R2.fq", fastqPipe)
	log.Printf("Finished Indexing Genome...\n")
	start := time.Now()

	log.Printf("Starting alignment worker...\n")
	workerWaiter.Add(numWorkers)
	for i := 0; i < numWorkers; i++ {
		go RoutineSimplyGiraf(&genome, tiles, tileSize, stepSize, scoreMatrix, fastqPipe, girafPipe, &workerWaiter)
	}

	go giraf.GirafPairChanToFile(output, girafPipe, &writerWaiter)
	writerWaiter.Add(1)
	workerWaiter.Wait()
	close(girafPipe)
	log.Printf("Aligners finished and channel closed\n")
	writerWaiter.Wait()
	log.Printf("Sam writer finished and we are all done\n")
	stop := time.Now()
	duration := stop.Sub(start)
	log.Printf("Aligned %d reads in %s (%.1f reads per second).\n", len(simReads)*2, duration, float64(len(simReads)*2)/duration.Seconds())
}*/

func BenchmarkTestWorkerWithWriting(b *testing.B) {
	if *cpuprofile != "" {
		f, err := os.Create(*cpuprofile)
		if err != nil {
			log.Fatal("could not create CPU profile: ", err)
		}
		defer f.Close() // error handling omitted for example
		if err := pprof.StartCPUProfile(f); err != nil {
			log.Fatal("could not start CPU profile: ", err)
		}
		defer pprof.StopCPUProfile()
	}
	//for n := 0; n < b.N; n++ {
	b.ResetTimer()
	var output string = "/dev/stdout"
	var tileSize int = 32
	var stepSize int = 32
	var numberOfReads int = 5000
	var readLength int = 150
	var mutations int = 2
	var workerWaiter, writerWaiter sync.WaitGroup
	var numWorkers int = 8
	var scoreMatrix = HumanChimpTwoScoreMatrix
	genome := SimplyRead("testdata/rabsTHREEspine.fa")

	log.Printf("Simulating reads...\n")
	simReads := RandomPairedReads(genome, readLength, numberOfReads, mutations)
	os.Remove("testdata/simReads_R1.fq")
	os.Remove("testdata/simReads_R2.fq")
	fastq.WritePair("testdata/simReads_R1.fq", "testdata/simReads_R2.fq", simReads)

	log.Printf("Reading in the genome (simple graph)...\n")
	log.Printf("Indexing the genome...\n")
	log.Printf("Making fastq channel...\n")
	fastqPipe := make(chan *fastq.PairedEndBig, 2408)

	log.Printf("Making giraf channel...\n")
	girafPipe := make(chan *giraf.GirafPair, 2408)

	tiles := IndexGenomeIntoMap(genome.Nodes, tileSize, stepSize)

	go readFastqGsw("testdata/simReads_R1.fq", "testdata/simReads_R2.fq", fastqPipe)
	log.Printf("Finished Indexing Genome...\n")
	start := time.Now()

	log.Printf("Starting alignment worker...\n")
	workerWaiter.Add(numWorkers)
	//matrixPool := NewMatrixPool()
	for i := 0; i < numWorkers; i++ {
		go RoutineSimplyGiraf(genome, tiles, tileSize, stepSize, scoreMatrix, fastqPipe, girafPipe, &workerWaiter)
	}

	go WriteSimpleGirafPair(output, girafPipe, &writerWaiter)
	writerWaiter.Add(1)
	workerWaiter.Wait()
	close(girafPipe)
	log.Printf("Aligners finished and channel closed\n")
	writerWaiter.Wait()
	log.Printf("Sam writer finished and we are all done\n")
	stop := time.Now()
	duration := stop.Sub(start)
	log.Printf("Aligned %d reads in %s (%.1f reads per second).\n", len(simReads)*2, duration, float64(len(simReads)*2)/duration.Seconds())
	//}
	if *memprofile != "" {
		f, err := os.Create(*memprofile)
		if err != nil {
			log.Fatal("could not create memory profile: ", err)
		}
		defer f.Close() // error handling omitted for example
		runtime.GC()    // get up-to-date statistics
		if err := pprof.WriteHeapProfile(f); err != nil {
			log.Fatal("could not write memory profile: ", err)
		}
	}
	
}

/*
func TestSimpleio(t *testing.T) {
	//var output string = "testdata/simple.giraf"
	var tileSize int = 32
	var stepSize int = 32
	var numberOfReads int = 1000
	var readLength int = 150
	var mutations int = 1
	var workerWaiter, writerWaiter sync.WaitGroup
	var numWorkers int = 8
	var scoreMatrix = HumanChimpTwoScoreMatrix
	genome := Read("testdata/gasAcu1.fa")

	log.Printf("Simulating reads...\n")
	simReads := RandomPairedReads(genome, readLength, numberOfReads, mutations)
	os.Remove("testdata/simReads_R1.fq")
	os.Remove("testdata/simReads_R2.fq")
	fastq.WritePair("testdata/simReads_R1.fq", "testdata/simReads_R2.fq", simReads)

	log.Printf("Reading in the genome (simple graph)...\n")
	log.Printf("Indexing the genome...\n")
	log.Printf("Making fastq channel...\n")
	fastqPipe := make(chan *fastq.PairedEndBig, 2408)

	log.Printf("Making giraf channel...\n")
	girafPipe := make(chan *giraf.GirafPair, 2408)

	tiles := IndexGenomeIntoMap(genome.Nodes, tileSize, stepSize)
	go fastq.ReadPairBigToChan("testdata/simReads_R1.fq", "testdata/simReads_R2.fq", fastqPipe)
	log.Printf("Finished Indexing Genome...\n")
	start := time.Now()

	log.Printf("Starting alignment worker...\n")
	workerWaiter.Add(numWorkers)
	for i := 0; i < numWorkers; i++ {
		go RoutineSimplyGiraf(genome, tiles, tileSize, stepSize, scoreMatrix, fastqPipe, girafPipe, &workerWaiter)
	}
	//go giraf.GirafPairChanToFile(output, girafPipe, &writerWaiter)
	go isGirafPairCorrect(girafPipe, genome, &writerWaiter, numberOfReads*2)
	writerWaiter.Add(1)
	workerWaiter.Wait()
	close(girafPipe)
	log.Printf("Aligners finished and channel closed\n")
	writerWaiter.Wait()
	//b.StopTimer()
	log.Printf("Giraf writer finished and we are all done\n")
	stop := time.Now()
	duration := stop.Sub(start)
	log.Printf("Aligned %d reads in %s (%.1f reads per second).\n", len(simReads)*2, duration, float64(len(simReads)*2)/duration.Seconds())
}*/
/*
func BenchmarkSimpleio(b *testing.B) {
	b.ReportAllocs()
	b.ResetTimer()
	var output string = "testdata/pairedTest.giraf"
	var tileSize int = 32
	var stepSize int = 32
	var numberOfReads int = 20000
	var readLength int = 150
	var mutations int = 1
	var workerWaiter, writerWaiter sync.WaitGroup
	var numWorkers int = 8
	var scoreMatrix = HumanChimpTwoScoreMatrix
	genome := ReadGenomeGraph("testdata/gasAcu1.fa")

	log.Printf("Simulating reads...\n")
	simReads := RandomPairedReads(genome, readLength, numberOfReads, mutations)
	os.Remove("testdata/simReads_R1.fq")
	os.Remove("testdata/simReads_R2.fq")
	fastq.WritePair("testdata/simReads_R1.fq", "testdata/simReads_R2.fq", simReads)

	if *cpuprofile != "" {
		f, err := os.Create(*cpuprofile)
		if err != nil {
			log.Fatal("could not create CPU profile: ", err)
		}
		defer f.Close() // error handling omitted for example
		if err := pprof.StartCPUProfile(f); err != nil {
			log.Fatal("could not start CPU profile: ", err)
		}
		defer pprof.StopCPUProfile()
	}
	log.Printf("Reading in the genome (simple graph)...\n")
	log.Printf("Indexing the genome...\n")
	log.Printf("Making fastq channel...\n")
	fastqPipe := make(chan *fastq.PairedEndBig, 2408)

	log.Printf("Making sam channel...\n")
	girafPipe := make(chan *giraf.GirafPair, 2408)

	tiles := IndexGenomeIntoMap(genome.Nodes, tileSize, stepSize)
	go fastq.ReadPairBigToChan("testdata/simReads_R1.fq", "testdata/simReads_R2.fq", fastqPipe)
	log.Printf("Finished Indexing Genome...\n")
	start := time.Now()

	log.Printf("Starting alignment worker...\n")
	workerWaiter.Add(numWorkers)
	for i := 0; i < numWorkers; i++ {
		go RoutineSimplyGiraf(genome, tiles, tileSize, stepSize, scoreMatrix, fastqPipe, girafPipe, &workerWaiter)
	}
	go giraf.GirafPairChanToFile(output, girafPipe, &writerWaiter)
	writerWaiter.Add(1)
	workerWaiter.Wait()
	close(girafPipe)
	log.Printf("Aligners finished and channel closed\n")
	writerWaiter.Wait()
	//b.StopTimer()
	log.Printf("Sam writer finished and we are all done\n")
	stop := time.Now()
	duration := stop.Sub(start)
	log.Printf("Aligned %d reads in %s (%.1f reads per second).\n", len(simReads)*2, duration, float64(len(simReads)*2)/duration.Seconds())

	if *memprofile != "" {
		f, err := os.Create(*memprofile)
		if err != nil {
			log.Fatal("could not create memory profile: ", err)
		}
		defer f.Close() // error handling omitted for example
		runtime.GC()    // get up-to-date statistics
		if err := pprof.WriteHeapProfile(f); err != nil {
			log.Fatal("could not write memory profile: ", err)
		}
	}
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
	girafPipe := make(chan *giraf.Giraf, 1)

	log.Printf("Starting alignment worker...\n")
	go RoutineFqToGiraf(genome, tiles, scoreMatrix, tileSize, stepSize, fastqPipe, girafPipe, &dummyWaiter)

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
