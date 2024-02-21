package genomeGraph

import (
	"flag"
	"log"
	"os"
	"runtime"
	"runtime/pprof"
	"strings"
	"sync"
	"testing"
	"time"

	"github.com/vertgenlab/gonomics/cigar"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/fastq"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/giraf"
	"github.com/vertgenlab/gonomics/numbers/parse"
)

var cpuprofile = flag.String("cpuprofile", "cpuprofile", "write cpu profile to `file`")
var memprofile = flag.String("memprofile", "memprofile", "write memory profile to `file`")

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
	b.ReportAllocs()
	var output string = "testdata/pairedTest.giraf"
	var tileSize int = 32
	var stepSize int = 32
	var numberOfReads int = 500
	var readLength int = 150
	var mutations int = 0
	var workerWaiter, writerWaiter sync.WaitGroup
	var numWorkers int = 8
	var scoreMatrix = HumanChimpTwoScoreMatrix
	genome := Read("testdata/bigGenome.sg")
	log.Printf("Reading in the genome (simple graph)...\n")
	log.Printf("Indexing the genome...\n")
	log.Printf("Making fastq channel...\n")
	fastqPipe := make(chan fastq.PairedEndBig, 824)

	log.Printf("Making sam channel...\n")
	samPipe := make(chan giraf.GirafPair, 824)

	log.Printf("Simulating reads...\n")
	simReads := RandomPairedReads(genome, readLength, numberOfReads, mutations)
	//os.Remove("testdata/simReads_R1.fq")
	//os.Remove("testdata/simReads_R2.fq")
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
	fileio.EasyRemove("testdata/simReads_R1.fq")
	fileio.EasyRemove("testdata/simReads_R2.fq")
	fileio.EasyRemove("testdata/pairedTest.giraf")

	//log.Printf("Aligned %d reads in %s (%.1f reads per second).\n", len(simReads)*2, duration, float64(len(simReads)*2)/duration.Seconds())
	log.Printf("Aligned %d reads in %s (%.1f reads per second).\n", numberOfReads*2, duration, float64(numberOfReads*2)/duration.Seconds())
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

func checkAlignment(aln giraf.Giraf, genome *GenomeGraph) bool {
	qName := strings.Split(aln.QName, "_")
	//if len(qName) < 5 {
	//	log.Fatalf("Error: input giraf file does not match simulation format...\n")
	//}
	if len(aln.Cigar) < 1 {
		return false
	}

	targetStart := aln.Path.TStart
	targetEnd := aln.Path.TEnd
	//if len(aln.Aln) < 1 {
	if aln.Cigar[0].Op == 'S' {
		//log.Printf("%s\n", giraf.GirafToString(aln))
		targetStart = targetStart - int(aln.Cigar[0].RunLen)
	}
	if aln.Cigar[len(aln.Cigar)-1].Op == 'S' {
		targetEnd = targetEnd + int(aln.Cigar[len(aln.Cigar)-1].RunLen)

		//}
	}
	if parse.StringToInt(qName[0]) == int(aln.Path.Nodes[0]) && parse.StringToInt(qName[1]) == targetStart && targetEnd == parse.StringToInt(qName[3]) {
		//log.Printf("%s\n", giraf.GirafToString(aln))
		//log.Printf("Results: %d != %d or %d != %d\n", headNode, aln.Path.Nodes[0], startPos, aln.Path.TStart)
		//	log.Printf("%s\n", giraf.GirafToString(aln))
		return true
	} else {
		//log.Printf("endPos=%d, right side cigar runLength: %d\n", endPos, aln.Aln[len(aln.Aln)-1].RunLen)
		//log.Printf("%s\n", giraf.GirafToString(aln))
		//log.Printf("Error: this read is not aligning correctly...\n")
	}
	return false
}
func percentOfFloat(part int, total int) float64 {
	return (float64(part) * float64(100)) / float64(total)
}

func isGirafPairCorrect(input <-chan giraf.GirafPair, genome *GenomeGraph, wg *sync.WaitGroup, numReads int) {
	var unmapped int = 0
	for pair := range input {
		if !checkAlignment(pair.Fwd, genome) {
			unmapped++
			//log.Printf("Error: failed alignment simulation...\n")
			//buf := GirafStringBuilder(pair.Fwd,&bytes.Buffer{})
			//log.Printf("%s\n", buf.String())
			//log.Printf("%s\n", giraf.GirafToString(pair.Rev))
		}

		if !checkAlignment(pair.Rev, genome) {
			//log.Printf("Error: failed alignment simulation...\n")
			//buf := GirafStringBuilder(pair.Fwd,&bytes.Buffer{})
			//log.Printf("%s\n", buf.String())
			unmapped++
		}
	}

	log.Printf("Mapped %d out of %d\n", numReads-unmapped, numReads)
	log.Printf("%f of the reads are mapping correctly\n", percentOfFloat(numReads-unmapped, numReads))

	wg.Done()
}

func TestAddSClip(t *testing.T) {
	// Test case 1: runLen < lengthOfRead
	cig := []cigar.Cigar{
		{RunLength: 10, Op: 'M'},
		{RunLength: 5, Op: 'D'},
		{RunLength: 15, Op: 'M'},
	}
	expected := []cigar.Cigar{
		{RunLength: 5, Op: 'S'},
		{RunLength: 10, Op: 'M'},
		{RunLength: 5, Op: 'D'},
		{RunLength: 15, Op: 'M'},
	}
	result := AddSClip(5, 30, cig)
	if !cigarSliceEqual(result, expected) {
		t.Errorf("Test case 1 failed: expected %v, got %v", expected, result)
	}

	// Test case 2: runLen >= lengthOfRead
	cig = []cigar.Cigar{
		{RunLength: 10, Op: 'M'},
		{RunLength: 5, Op: 'D'},
		{RunLength: 15, Op: 'M'},
	}
	expected = cig
	result = AddSClip(5, 20, cig)
	if !cigarSliceEqual(result, expected) {
		t.Errorf("Test case 2 failed: expected %v, got %v", expected, result)
	}
}

func cigarSliceEqual(a, b []cigar.Cigar) bool {
	if len(a) != len(b) {
		return false
	}
	for i := range a {
		if a[i] != b[i] {
			return false
		}
	}
	return true
}

func TestPerfectMatchBig(t *testing.T) {

	// Test case 1: Perfect match
	read := fastq.FastqBig{
		Seq: dna.StringToBases("ACGT"),
	}
	expected := int64(100)
	result := perfectMatchBig(read, HumanChimpTwoScoreMatrix)
	if result != expected {
		t.Errorf("Test case 1 failed: expected %d, got %d", expected, result)
	}

	// Test case 2: Reverse complement match
	read = fastq.FastqBig{
		Seq: dna.StringToBases("TGCA"),
	}
	expected = int64(100)
	result = perfectMatchBig(read, HumanChimpTwoScoreMatrix)
	if result != expected {
		t.Errorf("Test case 2 failed: expected %d, got %d", expected, result)
	}

	// Test case 3: Mismatch
	read = fastq.FastqBig{
		Seq: dna.StringToBases("ACTT"),
	}

	expected = int64(100)
	result = perfectMatchBig(read, HumanChimpTwoScoreMatrix)
	if result != expected {
		t.Errorf("Test case 3 failed: expected %d, got %d", expected, result)
	}

	// Test case 3: Mismatch
	read = fastq.FastqBig{
		Seq: dna.StringToBases("AAA"),
	}

	expected = int64(90)
	result = perfectMatchBig(read, HumanChimpTwoScoreMatrix)
	if result != expected {
		t.Errorf("Test case 4 failed: expected %d, got %d", expected, result)
	}
}
