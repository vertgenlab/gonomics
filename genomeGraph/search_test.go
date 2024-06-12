package genomeGraph

import (
	"strings"
	"testing"
	"time"

	"github.com/vertgenlab/gonomics/align"
	"github.com/vertgenlab/gonomics/cigar"
	"github.com/vertgenlab/gonomics/fastq"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/giraf"
	"github.com/vertgenlab/gonomics/numbers/parse"
)

func BenchmarkGirafAlignment(b *testing.B) {
	b.ResetTimer()
	b.ReportAllocs()

	var numberOfReads int = 10
	var readLength int = 150
	var mutations int = 1

	config := &GraphSettings{
		ScoreMatrix:    align.HumanChimpTwoScoreMatrix,
		GapPenalty:     -600,
		OpenGapPenalty: -150,
		TileSize:       32,
		StepSize:       32,
	}

	genome := Read("testdata/bigGenome.sg")
	tiles := IndexGenomeIntoMap(genome.Nodes, config.TileSize, config.StepSize)

	simReads := RandomPairedReads(genome, readLength, numberOfReads, mutations)
	fqOne := "testdata/simReads_R1.fq"
	fqTwo := "testdata/simReads_R2.fq"

	fastq.WritePair(fqOne, fqTwo, simReads)
	fqs := []fastq.FastqBig{}
	file := fileio.NewByteReader(fqOne)

	for read, done := fastq.ReadFqBig(file); !done; read, done = fastq.ReadFqBig(file) {
		fqs = append(fqs, read)
	}

	matrix := NewSwMatrix(defaultMatrixSize)
	seedPool := NewMemSeedPool()
	dnaPool := NewDnaPool()
	seedBuildHelper := newSeedBuilder()
	scorekeeper := scoreKeeper{}
	dynamicKeeper := dynamicScoreKeeper{}
	var correct int
	start := time.Now()
	for n := 0; n < b.N; n++ {
		var results []giraf.Giraf = make([]giraf.Giraf, len(fqs))
		correct = 0
		for i := 0; i < len(fqs); i++ {
			results[i] = *GraphSmithWatermanToGiraf(genome, fqs[i], tiles, config, &matrix, &seedPool, &dnaPool, scorekeeper, dynamicKeeper, seedBuildHelper)
			if isCorrectCoord(&results[i]) {
				correct++
			}
		}

	}
	stop := time.Now()
	duration := stop.Sub(start)
	b.Logf("Mapping results: %d / %d reads correctly", correct, len(fqs))
	b.Logf("Aligned %d reads in %s (%.1f reads per second).\n", len(simReads), duration, float64(len(simReads)*2)/duration.Seconds())
	fileio.EasyRemove("testdata/simReads_R1.fq")
	fileio.EasyRemove("testdata/simReads_R2.fq")
}

func isCorrectCoord(result *giraf.Giraf) bool {
	name := strings.Split(result.QName, "_")
	// TODO: Need better logic to look at node
	return parse.StringToInt(name[1]) == result.Path.TStart-result.QStart && parse.StringToInt(name[3]) == result.Path.TEnd && cigar.QueryRunLen(result.Cigar) == 150
}
