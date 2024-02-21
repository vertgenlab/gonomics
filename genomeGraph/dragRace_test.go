package genomeGraph

import (
	"os"
	"runtime/pprof"
	"sync"
	"testing"
	"time"

	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fastq"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/sam"
)

const Profile int = 0

func TestQuickMemPool(t *testing.T) {
	var tileSize int = 14
	var stepSize int = 4
	var numberOfReads int = 20
	var readLength int = 150
	var mutations int = 0
	var workerWaiter, writerWaiter sync.WaitGroup
	var numWorkers int = 1
	var scoreMatrix = HumanChimpTwoScoreMatrix

	t.Logf("Reading in the genome (simple graph)...\n")

	genome := Read("testdata/mini.gg")

	t.Logf("Indexing the genome...\n")
	tiles := IndexGenomeIntoMap(genome.Nodes, tileSize, stepSize)

	fastqPipe := make(chan fastq.FastqBig, 824)

	samPipe := make(chan sam.Sam, 824)

	simReads := RandomReads(genome, readLength, numberOfReads, mutations)
	fastq.Write("testdata/simReads.fq", simReads)

	go fastq.ReadBigToChan("testdata/simReads.fq", fastqPipe)
	writerWaiter.Add(1)
	go sam.WriteFromChan(samPipe, "testdata/sim.sam", sam.Header{}, &writerWaiter)

	time.Sleep(5 * time.Second)

	if Profile > 0 {
		f, err := os.Create("testdata/cpuprofile.data")
		exception.PanicOnErr(err)
		defer f.Close()
		err = pprof.StartCPUProfile(f)
		exception.PanicOnErr(err)
	}
	start := time.Now()
	workerWaiter.Add(numWorkers)
	for i := 0; i < numWorkers; i++ {
		go gswWorkerMemPool(genome, tiles, tileSize, stepSize, scoreMatrix, fastqPipe, samPipe, &workerWaiter)
	}

	workerWaiter.Wait()
	stop := time.Now()
	pprof.StopCPUProfile()
	close(samPipe)

	writerWaiter.Wait()

	duration := stop.Sub(start)
	t.Logf("Aligned %d reads in %s (%.1f reads per second).\n", len(simReads), duration, float64(len(simReads))/duration.Seconds())

	fileio.EasyRemove("testdata/simReads.fq")
	fileio.EasyRemove("testdata/sim.sam")
}
