package genomeGraph

import (
	"github.com/vertgenlab/gonomics/common"
	"github.com/vertgenlab/gonomics/fastq"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/sam"
	"log"
	"os"
	"runtime/pprof"
	"sync"
	"testing"
	"time"
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

	log.Printf("Reading in the genome (simple graph)...\n")
	//genome := Read("testdata/bigGenome.sg")
	//genome := Read("testdata/rabsBepaChrI.gg")
	//genome := Read("testdata/tiny.gg")
	genome := Read("testdata/mini.gg")

	log.Printf("Indexing the genome...\n")
	tiles := IndexGenomeIntoMap(genome.Nodes, tileSize, stepSize)

	log.Printf("Making fastq channel...\n")
	fastqPipe := make(chan *fastq.FastqBig, 824)

	log.Printf("Making sam channel...\n")
	samPipe := make(chan *sam.Aln, 824)

	log.Printf("Simulating reads...\n")
	simReads := RandomReads(genome, readLength, numberOfReads, mutations)
	fastq.Write("testdata/simReads.fq", simReads)

	//header := NodesHeader(genome.Nodes)

	go fastq.ReadBigToChan("testdata/simReads.fq", fastqPipe)
	writerWaiter.Add(1)
	go sam.SamChanToFile(samPipe, "testdata/sim.sam", nil, &writerWaiter)

	log.Printf("Starting alignment worker...\n")
	time.Sleep(5 * time.Second)

	if Profile > 0 {
		f, err := os.Create("testdata/cpuprofile.data")
		common.ExitIfError(err)
		defer f.Close()
		err = pprof.StartCPUProfile(f)
		common.ExitIfError(err)
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

	log.Printf("Aligners finished and channel closed\n")
	writerWaiter.Wait()
	log.Printf("Sam writer finished and we are all done\n")
	duration := stop.Sub(start)
	log.Printf("Aligned %d reads in %s (%.1f reads per second).\n", len(simReads), duration, float64(len(simReads))/duration.Seconds())
	fileio.EasyRemove("testdata/simReads.fq")
	fileio.EasyRemove("testdata/sim.sam")
}
