package genomeGraph

import (
	"os"
	"testing"
	"time"

	"github.com/vertgenlab/gonomics/align"
	"github.com/vertgenlab/gonomics/fastq"
	"github.com/vertgenlab/gonomics/fileio"
)

func TestGraphSmithWatermanToGiraf(t *testing.T) {
	// var output string = "testdata/pairedTest.giraf"
	var tileSize int = 32
	var stepSize int = 32
	var scoreMatrix = align.HumanChimpTwoScoreMatrix
	var numberOfReads int = 1000
	var mutations int = 2

	config := &GraphSettings{
		ScoreMatrix: scoreMatrix,
		GapPenalty:  -600,
		TileSize:    tileSize,
		StepSize:    stepSize,
	}
	memoryPool := MatrixPoolMemory(defaultMatrixSize)

	genome := Read("testdata/bigGenome.sg")
	simReads := RandomPairedReads(genome, 150, numberOfReads, mutations)
	fqOneSimFile, fqTwoSimFile := "testdata/simulated-reads_R1.fq", "testdata/simulated-reads_R2.fq"
	fastq.WritePair(fqOneSimFile, fqTwoSimFile, simReads)

	tiles := IndexGenomeIntoMap(genome.Nodes, tileSize, stepSize)

	fileOne := fileio.NewByteReader(fqOneSimFile)
	fqs := []fastq.FastqBig{}
	var counter int = 0

	seedPool := NewMemSeedPool()

	seedBuildHelper := newSeedBuilder()
	scorekeeper := scoreKeeper{}

	for read, done := fastq.ReadFqBig(fileOne); !done; read, done = fastq.ReadFqBig(fileOne) {
		fqs = append(fqs, read)
	}
	start := time.Now()
	for i := 0; i < len(fqs); i++ {
		result := GraphSmithWatermanToGiraf(genome, fqs[i], tiles, config, memoryPool, &seedPool, scorekeeper, seedBuildHelper)
		if IsCorrectCoord(result) {
			counter++
		}
	}
	stop := time.Now()
	duration := stop.Sub(start)
	t.Logf("Aligned %d reads in %s (%.1f reads per second).\n", len(fqs), duration, float64(len(fqs))/duration.Seconds())
	t.Logf("Aligned %d/%d reads correctly!\n", counter, len(fqs))
	os.Remove(fqOneSimFile)
	os.Remove(fqTwoSimFile)
}
