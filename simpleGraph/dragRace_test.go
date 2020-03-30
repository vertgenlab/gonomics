package simpleGraph

import (
	"github.com/vertgenlab/gonomics/common"
	"github.com/vertgenlab/gonomics/fastq"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/sam"
	"log"
	"os"
	"runtime/pprof"
	"strings"
	"sync"
	"testing"
	"time"
)

func TestQuickMemPool(t *testing.T) {
	var tileSize int = 32
	var stepSize int = 32
	var numberOfReads int = 10
	var readLength int = 150
	var mutations int = 0
	var workerWaiter, writerWaiter sync.WaitGroup
	var numWorkers int = 2
	var scoreMatrix = HumanChimpTwoScoreMatrix

	log.Printf("Reading in the genome (simple graph)...\n")
	genome, _ := Read("testdata/bigGenome.sg")
	//genome, _ := Read("testdata/rabsBepaChrI.gg")

	log.Printf("Indexing the genome...\n")
	tiles := indexGenomeIntoMap(genome.Nodes, tileSize, stepSize)

	log.Printf("Making fastq channel...\n")
	fastqPipe := make(chan *fastq.FastqBig, 824)

	log.Printf("Making sam channel...\n")
	samPipe := make(chan *sam.SamAln, 824)

	log.Printf("Simulating reads...\n")
	simReads := RandomReads(genome.Nodes, readLength, numberOfReads, mutations)
	fastq.Write("testdata/simReads.fq", simReads)

	header := NodesHeader(genome.Nodes)

	go fastq.ReadBigToChan("testdata/simReads.fq", fastqPipe)
	writerWaiter.Add(1)
	go sam.SamChanToFile(samPipe, "testdata/sim.sam", header, &writerWaiter)

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
	strictCheckFile("testdata/sim.sam", genome)
	log.Printf("Passed alignment check!!!\n")
}

func strictCheckAlignment(aln *sam.SamAln, genome *SimpleGraph) {
	qName := strings.Split(aln.QName, "_")
	if strings.Compare(genome.Nodes[common.StringToUint32(qName[0])].Name, aln.RName) == 0 && aln.Pos == common.StringToInt64(qName[1]) {
		// do nothing
	} else {
		log.Fatalf("Error: Incorrect Alignment\n")
	}
}

func strictCheckFile(filename string, genome *SimpleGraph) {
	var curr *sam.SamAln
	var done bool

	file := fileio.EasyOpen(filename)
	sam.ReadHeader(file)
	defer file.Close()

	for curr, done = sam.NextAlignment(file); done != true; curr, done = sam.NextAlignment(file) {
		strictCheckAlignment(curr, genome)
	}
}
