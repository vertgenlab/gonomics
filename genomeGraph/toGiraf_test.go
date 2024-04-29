package genomeGraph

// import (
// 	"log"
// 	"testing"
// 	"time"

// 	"github.com/vertgenlab/gonomics/align"
// 	"github.com/vertgenlab/gonomics/fastq"
// 	"github.com/vertgenlab/gonomics/fileio"
// 	"github.com/vertgenlab/gonomics/giraf"
// )

// func TestGraphSmithWatermanToGiraf(t *testing.T) {
// 	// var output string = "testdata/pairedTest.giraf"
// 	var tileSize int = 32
// 	var stepSize int = 32
// 	var scoreMatrix = align.HumanChimpTwoScoreMatrix
// 	config := &GraphSettings{
// 		ScoreMatrix: scoreMatrix,
// 		GapPenalty:  -600,
// 		TileSize:    tileSize,
// 		StepSize:    stepSize,
// 	}
// 	memoryPool := MatrixPoolMemory(defaultMatrixSize)

// 	genome := Read("testdata/bigGenome.sg")
// 	tiles := IndexGenomeIntoMap(genome.Nodes, tileSize, stepSize)

// 	fqOne := "testdata/simReads_R1.fq"

// 	fileOne := fileio.NewByteReader(fqOne)
// 	fqs := []fastq.FastqBig{}

// 	seedPool := NewMemSeedPool()

// 	seedBuildHelper := newSeedBuilder()
// 	scorekeeper := scoreKeeper{}

// 	for read, done := fastq.ReadFqBig(fileOne); !done; read, done = fastq.ReadFqBig(fileOne) {
// 		fqs = append(fqs, read)
// 	}
// 	start := time.Now()
// 	results := []giraf.Giraf{}
// 	for i := 0; i < len(fqs); i++ {
// 		results = append(results, *GraphSmithWatermanToGiraf(genome, fqs[i], tiles, config, memoryPool, &seedPool, scorekeeper, seedBuildHelper))
// 	}
// 	stop := time.Now()
// 	duration := stop.Sub(start)
// 	log.Printf("Aligned %d reads in %s (%.1f reads per second).\n", len(fqs), duration, float64(len(fqs))/duration.Seconds())
// 	var counter int = 0
// 	for j := 0; j < len(results); j++ {
// 		if IsCorrectCoord(results[j]) {
// 			counter++
// 		}
// 	}
// 	log.Printf("Aligned %d/%d reads correctly!\n", counter, len(results))
// }
