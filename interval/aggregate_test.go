package interval

import (
	"testing"

	"github.com/vertgenlab/gonomics/bed"
)

func TestMergeIntervals(t *testing.T) {
	var testIntervals []Interval = []Interval{
		&bed.Bed{Chrom: "1", ChromStart: 1, ChromEnd: 8},
		&bed.Bed{Chrom: "1", ChromStart: 3, ChromEnd: 5},
		&bed.Bed{Chrom: "1", ChromStart: 4, ChromEnd: 4},
		&bed.Bed{Chrom: "1", ChromStart: 5, ChromEnd: 9},
		&bed.Bed{Chrom: "1", ChromStart: 6, ChromEnd: 11},
		&bed.Bed{Chrom: "1", ChromStart: 7, ChromEnd: 7},
		&bed.Bed{Chrom: "1", ChromStart: 8, ChromEnd: 10},
		&bed.Bed{Chrom: "1", ChromStart: 9, ChromEnd: 10},
		&bed.Bed{Chrom: "1", ChromStart: 11, ChromEnd: 12},
		&bed.Bed{Chrom: "1", ChromStart: 13, ChromEnd: 14},
		&bed.Bed{Chrom: "2", ChromStart: 1, ChromEnd: 6},
	}
	answer := MergeIntervals(testIntervals)

	if len(answer) != 3 {
		t.Errorf("ERROR: Problem merging intervals")
	}
}
