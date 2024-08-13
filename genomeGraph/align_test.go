package genomeGraph

import (
	"strings"
	"sync"
	"testing"
	"time"

	"github.com/vertgenlab/gonomics/cigar"
	"github.com/vertgenlab/gonomics/fastq"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/giraf"
	"github.com/vertgenlab/gonomics/numbers/parse"
)

func BenchmarkGsw(b *testing.B) {
	b.ReportAllocs()
	var tileSize int = 32
	var stepSize int = 32
	var numberOfReads int = 100
	var readLength int = 150
	var mutations int = 1
	var workerWaiter, writerWaiter sync.WaitGroup
	var numWorkers int = 1

	config := DefaultAlignment()

	b.ResetTimer()
	genome := Read("testdata/bigGenome.sg")
	tiles := IndexGenomeIntoMap(genome.Nodes, tileSize, stepSize)

	simReads := RandomPairedReads(genome, readLength, numberOfReads, mutations)
	fqOne, fqTwo := "testdata/simReads_R1.fq", "testdata/simReads_R2.fq"
	fastq.WritePair(fqOne, fqTwo, simReads)
	//var output string = "pprof/pairedTest.giraf"

	fastqPipe, girafPipe := make(chan fastq.PairedEndBig, 2408), make(chan giraf.GirafPair, 2408)
	go readFastqGsw("testdata/simReads_R1.fq", "testdata/simReads_R2.fq", fastqPipe)

	for n := 0; n < b.N; n++ {
		start := time.Now()
		workerWaiter.Add(numWorkers)
		dp := NewMatrixPool(defaultMatrixSize)

		for i := 0; i < numWorkers; i++ {
			go RoutineFqPairToGiraf(genome, tiles, config, dp, fastqPipe, girafPipe, &workerWaiter)
		}
		//go giraf.GirafPairChanToFile(output, girafPipe, &writerWaiter)
		go isGirafPairCorrect(girafPipe, &writerWaiter, readLength, 2*len(simReads), b)
		writerWaiter.Add(1)

		stop := time.Now()
		duration := stop.Sub(start)
		//log.Printf("Aligned %d reads in %s (%.1f reads per second).\n", len(simReads)*2, duration, float64(len(simReads)*2)/duration.Seconds())
		b.Logf("Aligned %d reads in %s (%.1f reads per second).\n", 2*len(simReads), duration, float64(2*len(simReads))/duration.Seconds())
	}
	workerWaiter.Wait()
	close(girafPipe)
	writerWaiter.Wait()
	fileio.EasyRemove(fqOne)
	fileio.EasyRemove(fqTwo)
}

func checkAlignment(result giraf.Giraf, readLength int) bool {
	if len(result.Cigar) == 0 {
		return false
	}
	name := strings.Split(result.QName, "_")
	return parse.StringToInt(name[1]) == result.Path.TStart-result.QStart && parse.StringToInt(name[3]) == result.Path.TEnd+(readLength-result.QEnd-1-result.QStart) && cigar.QueryLength(result.Cigar) == readLength
}
func percentOfFloat(part int, total int) float64 {
	return (float64(part) * float64(100)) / float64(total)
}

func isGirafPairCorrect(input <-chan giraf.GirafPair, wg *sync.WaitGroup, readLength int, numReads int, b *testing.B) {
	var unmapped int = 0
	for pair := range input {
		//fmt.Printf("%s\n", cigar.ToString(pair.Fwd.Cigar))
		//fmt.Printf("%s\n%s\n", cigar.ToString(pair.Fwd.Cigar), giraf.GirafToString(&pair.Fwd))
		if !checkAlignment(pair.Fwd, readLength) {
			unmapped++
			// if !cigar.IsUnmapped(pair.Fwd.Cigar) {
			// 	if !cigar.AllEqual(pair.Fwd.Cigar, []cigar.Cigar{{150, cigar.Match}}) {
			// 		fmt.Printf("%s\n%s\n", cigar.ToString(pair.Fwd.Cigar), giraf.GirafToString(&pair.Fwd))
			// 	}
			// }

			// name := strings.Split(pair.Fwd.QName, "_")
			// if len(pair.Fwd.Cigar)> 0&& parse.StringToInt(name[3]) != pair.Fwd.Path.TEnd +(readLength-pair.Fwd.QEnd-1-pair.Fwd.QStart) {
			// 	b.Errorf("read==%s\n", giraf.GirafToString(&pair.Fwd))
			// 	// b.Errorf("%d != %d\n", parse.StringToInt(name[3]), pair.Fwd.Path.TEnd + (readLength-pair.Fwd.QEnd-1-pair.Fwd.QStart) )
			// }

		}
		if !checkAlignment(pair.Rev, readLength) {
			unmapped++
		}
	}
	b.Logf("Mapped %d out of %d\n", numReads-unmapped, numReads)
	b.Logf("%f of the reads are mapping correctly\n", percentOfFloat(numReads-unmapped, numReads))
	wg.Done()
}
