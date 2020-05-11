package simpleGraph
import (
	"github.com/vertgenlab/gonomics/fastq"
	"github.com/vertgenlab/gonomics/giraf"
	"testing"
	"log"
	"sync"

	"os"
)

func TestViewGirafAlign(t *testing.T) {
	var tileSize int = 32
	var stepSize int = 32
	var numberOfReads int = 100
	var readLength int = 150
	var mutations int = 0
	var workerWaiter sync.WaitGroup
	//, writerWaiter
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
	//start := time.Now()

	log.Printf("Starting alignment worker...\n")
	workerWaiter.Add(numWorkers)
	for i := 0; i < numWorkers; i++ {
		go RoutineFqPairToGiraf(genome, tiles, tileSize, stepSize, scoreMatrix, fastqPipe, samPipe, &workerWaiter)
	}
	for alignPair := range samPipe {
		log.Printf("forward:\n%s\n", ViewGirafAlign(alignPair.Fwd, genome))
		log.Printf("reverse:\n%s\n", ViewGirafAlign(alignPair.Rev, genome))
	}
}
/*
import (
	"github.com/vertgenlab/gonomics/chromInfo"
	"github.com/vertgenlab/gonomics/fastq"
	"github.com/vertgenlab/gonomics/giraf"
	"github.com/vertgenlab/gonomics/sam"
	//"github.com/vertgenlab/gonomics/align"
	"log"
	"os"
	"sync"
	"testing"
	"time"
)


//TODO: finish writing matrix helper function
func TestScoreMatrixHelper(t *testing.T) {
	var scoreMatrixSlice [][][]int64 = [][][]int64{align.HumanChimpTwoScoreMatrix, align.HoxD55ScoreMatrix, align.MouseRatScoreMatrix}
	for i := 0; i < len(scoreMatrixSlice);i++ {
		help := getScoreMatrixHelp(scoreMatrixSlice[i])
		maxMatch, minMatch, leastSevereMismatch, leastSevereMatchMismatchChange := help.MaxMatch, help.MinMatch, help.LeastSevereMismatch, help.LeastSevereMatchMismatchChange
		log.Printf("maxMatch=%d, minMatch=%d, leastSevereMismatch=%d, leastSevereMatchMismatchChange=%d", maxMatch, minMatch, leastSevereMismatch, leastSevereMatchMismatchChange)
	}
}*/
