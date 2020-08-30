package simpleGraph

import (
	"flag"
	"github.com/vertgenlab/gonomics/common"
	"github.com/vertgenlab/gonomics/fastq"
	"github.com/vertgenlab/gonomics/giraf"
	"log"
	"os"
	"runtime"
	"runtime/pprof"
	"strings"
	"sync"
	"testing"
	"time"
)

var cpuprofile = flag.String("cpuprofile", "", "write cpu profile to `file`")
var memprofile = flag.String("memprofile", "", "write memory profile to `file`")

func BenchmarkGsw(b *testing.B) {

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
	gorout, err := os.Create("/dev/stdout")
	if err != nil {
		log.Fatal("could not create GOroutines profile: ", err)
	}
	defer gorout.Close()

	b.ReportAllocs()
	//var output string = "/dev/stdout"
	//var output string = "testdata/rabs_test.giraf"
	var tileSize int = 32
	var stepSize int = 32
	var numberOfReads int = 25000
	var readLength int = 150
	var mutations int = 1
	var workerWaiter, writerWaiter sync.WaitGroup
	var numWorkers int = 4
	var scoreMatrix = HumanChimpTwoScoreMatrix
	b.ResetTimer()
	genome := Read("testdata/rabsTHREEspine.fa")
	log.Printf("Reading in the genome (simple graph)...\n")
	log.Printf("Indexing the genome...\n")

	fastqPipe := make(chan fastq.PairedEndBig, 2408)
	girafPipe := make(chan giraf.GirafPair, 2408)

	log.Printf("Simulating reads...\n")
	simReads := RandomPairedReads(genome, readLength, numberOfReads, mutations)
	os.Remove("testdata/simReads_R1.fq")
	os.Remove("testdata/simReads_R2.fq")

	fastq.WritePair("testdata/simReads_R1.fq", "testdata/simReads_R2.fq", simReads)

	tiles := IndexGenomeIntoMap(genome.Nodes, tileSize, stepSize)

	go readFastqGsw("testdata/simReads_R1.fq", "testdata/simReads_R2.fq", fastqPipe)
	log.Printf("Finished Indexing Genome...\n")

	log.Printf("Starting alignment worker...\n")
	workerWaiter.Add(numWorkers)
	start := time.Now()
	for i := 0; i < numWorkers; i++ {
		go RoutineFqPairToGiraf(genome, tiles, tileSize, stepSize, scoreMatrix, fastqPipe, girafPipe, &workerWaiter)
	}
	//go SimpleWriteGirafPair(output, girafPipe, &writerWaiter)
	go isGirafPairCorrect(girafPipe, genome, &writerWaiter, 2*len(simReads))
	writerWaiter.Add(1)
	workerWaiter.Wait()
	close(girafPipe)

	writerWaiter.Wait()

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

	//log.Printf("Aligned %d reads in %s (%.1f reads per second).\n", 50000*2, duration, float64(50000*2)/duration.Seconds())

}

func checkAlignment(aln giraf.Giraf, genome *SimpleGraph) bool {
	qName := strings.Split(aln.QName, "_")

	if len(aln.Cigar) < 1 {
		return false
	}

	targetStart := aln.Path.TStart
	targetEnd := aln.Path.TEnd

	if aln.Cigar[0].Op == 'S' {

		targetStart = targetStart - int(aln.Cigar[0].RunLen)
	}
	if aln.Cigar[len(aln.Cigar)-1].Op == 'S' {
		targetEnd = targetEnd + int(aln.Cigar[len(aln.Cigar)-1].RunLen)

	}
	if common.StringToInt(qName[0]) == int(aln.Path.Nodes[0]) && common.StringToInt(qName[1]) == targetStart && targetEnd == common.StringToInt(qName[3]) {

		return true
	} else {

	}
	return false
}
func percentOfFloat(part int, total int) float64 {
	return (float64(part) * float64(100)) / float64(total)
}

func isGirafPairCorrect(input <-chan giraf.GirafPair, genome *SimpleGraph, wg *sync.WaitGroup, numReads int) {
	var unmapped int = 0
	for pair := range input {
		if !checkAlignment(pair.Fwd, genome) {
			unmapped++
		}

		if !checkAlignment(pair.Rev, genome) {

			unmapped++
		}
	}
	log.Printf("Mapped %d out of %d\n", numReads-unmapped, numReads)
	log.Printf("%f of the reads are mapping correctly\n", percentOfFloat(numReads-unmapped, numReads))

	wg.Done()
}
