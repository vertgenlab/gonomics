package genomeGraph

import (
	"flag"
	"log"
	"os"
	"runtime/pprof"
	"strings"
	"sync"
	"testing"
	"time"

	"github.com/vertgenlab/gonomics/align"
	"github.com/vertgenlab/gonomics/fastq"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/giraf"
	"github.com/vertgenlab/gonomics/numbers/parse"
)

var cpuprofile = flag.String("cpuprofile", "cpuprofile", "write cpu profile to `file`")
var memprofile = flag.String("memprofile", "cpuprofile", "write memory profile to `file`")

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
	//var output string = "/dev/stdout"
	var output string = "testdata/rabs_test.giraf"
	var tileSize int = 32
	var stepSize int = 32
	var numberOfReads int = 1
	var readLength int = 150
	var mutations int = 1
	var workerWaiter, writerWaiter sync.WaitGroup
	var numWorkers int = 6
	settings := &GraphSettings{
		ScoreMatrix: align.HumanChimpTwoScoreMatrix,
		GapPenalty:  -600,
		TileSize:    tileSize,
		StepSize:    stepSize,
	}

	b.ResetTimer()
	genome := Read("testdata/bigGenome.sg")

	fastqPipe := make(chan fastq.PairedEndBig, 2408)
	girafPipe := make(chan giraf.GirafPair, 2408)
	simReads := RandomPairedReads(genome, readLength, numberOfReads, mutations)

	fastq.WritePair("testdata/simReads_R1.fq", "testdata/simReads_R2.fq", simReads)
	tiles := IndexGenomeIntoMap(genome.Nodes, tileSize, stepSize)

	go readFastqGsw("testdata/simReads_R1.fq", "testdata/simReads_R2.fq", fastqPipe)
	workerWaiter.Add(numWorkers)
	start := time.Now()
	for i := 0; i < numWorkers; i++ {
		go RoutineFqPairToGiraf(genome, tiles, settings.TileSize, settings.StepSize, settings.ScoreMatrix, fastqPipe, girafPipe, &workerWaiter)
	}
	go SimpleWriteGirafPair(output, girafPipe, &writerWaiter)
	writerWaiter.Add(1)
	workerWaiter.Wait()
	close(girafPipe)

	writerWaiter.Wait()

	stop := time.Now()
	duration := stop.Sub(start)
	//log.Printf("Aligned %d reads in %s (%.1f reads per second).\n", len(simReads)*2, duration, float64(len(simReads)*2)/duration.Seconds())
	log.Printf("Aligned %d reads in %s (%.1f reads per second).\n", len(simReads)*2, duration, float64(len(simReads)*2)/duration.Seconds())
	fileio.MustRemove(output)

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
