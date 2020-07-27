package main

import (
	"github.com/vertgenlab/gonomics/bed"
	"testing"
)

func TestSequelOverlap(t *testing.T) {
	init := &Settings{
		Input:           "testdata/test.vcf",
		Output:          "/dev/stdout",
		SelectFile:      "testdata/test.bed",
		Extend:          0,
		NonOverlap:      false,
		Threads:         1,
		PercentOverlap:  0,
		BaseOverlap:     0,
		Aggregate:       false,
		Relationship:    "any",
		MergedOutput:    false,
		SwapTargetQuery: false,
	}

	answer := sequelOverlap(init)

	for val := range answer {
		b := val.answer[0].(*bed.Bed)
		if b.Chrom != "chr1" || b.ChromStart != 100 || b.ChromEnd != 200 {
			t.Errorf("ERROR: Problem with sequelOverlap cmd")
		}
	}
}
