package genomeGraph

import (
	"testing"

	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/fastq"
)

func TestBlastSeed(t *testing.T) {
	seed := &SeedDev{
		TargetId:    0,
		TargetStart: 10,
		Length:      5,
	}

	read := fastq.Fastq{
		Name: "read1",
		Seq:  []dna.Base{dna.A, dna.T, dna.C, dna.G, dna.A},
	}

	scoreMatrix := [][]int64{
		{2, -1, -1, -1},
		{-1, 2, -1, -1},
		{-1, -1, 2, -1},
		{-1, -1, -1, 2},
	}

	result := BlastSeed(seed, read, scoreMatrix)
	expected := int64(10)

	if result != expected {
		t.Errorf("Error: Unexpected result, got: %d, expected: %d", result, expected)
	}
}
