package interval

import (
	"github.com/vertgenlab/gonomics/bed"
	"testing"
)

func TestSmallerIntervals(t *testing.T) {
	var testIntervals IntervalSlice = IntervalSlice{
		&bed.Bed{Chrom: "1", ChromStart: 2, ChromEnd: 4},
		&bed.Bed{Chrom: "1", ChromStart: 2, ChromEnd: 3},
		&bed.Bed{Chrom: "1", ChromStart: 0, ChromEnd: 1},
		&bed.Bed{Chrom: "1", ChromStart: 1, ChromEnd: 2},
		&bed.Bed{Chrom: "1", ChromStart: 4, ChromEnd: 4},
	}
	sortIntervals(testIntervals)
	newIntervals, smaller := computeSmallerIntervals(testIntervals)
	if smaller[newIntervals[2]][0].GetChromStart() != 2 ||
		smaller[newIntervals[2]][0].GetChromEnd() != 3 {
		t.Errorf("ERROR: Problem computing smaller intervals")
	}
}
