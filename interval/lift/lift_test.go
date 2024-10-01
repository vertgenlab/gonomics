package lift

import (
	"testing"

	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/interval"
)

var ReadingTests = []struct {
	Infile   string
	Expected []Lift
}{
	{"testdata/elements1.bed", []Lift{bed.Bed{Chrom: "chr1", ChromStart: 2, ChromEnd: 10}, bed.Bed{Chrom: "chr1", ChromStart: 13, ChromEnd: 15}, bed.Bed{Chrom: "chr1", ChromStart: 50, ChromEnd: 500}}},
}

func TestReading(t *testing.T) {
	var l []Lift
	for _, v := range ReadingTests {
		l = GoRead(v.Infile)
		if !allAreEqual(l, v.Expected) {
			t.Errorf("Error in lift GoRead. Input did not match expected.")
		}
	}
}

func TestIntervalSliceToLift(t *testing.T) {
	// Create intervals where some implement Lift, others don't
	intervals := []interval.Interval{
		bed.Bed{Chrom: "chr1", ChromStart: 1, ChromEnd: 30},
		bed.Bed{Chrom: "chr2", ChromStart: 300, ChromEnd: 400},
	}
	lifts := IntervalSliceToLift(intervals)
	expected := []Lift{
		bed.Bed{Chrom: "chr1", ChromStart: 1, ChromEnd: 30},
		bed.Bed{Chrom: "chr2", ChromStart: 300, ChromEnd: 400},
	}
	if !allAreEqual(lifts, expected) {
		t.Errorf("Error: IntervalSliceToLift(intervals, expected) do not equal.\n")
	}
}

func TestMatchOverlapLen(t *testing.T) {
	var one bed.Bed = bed.Bed{Chrom: "chr1", ChromStart: 5, ChromEnd: 15}
	var two bed.Bed = bed.Bed{Chrom: "chr1", ChromStart: 11, ChromEnd: 14}
	overlapSize := MatchOverlapLen(one.GetChromStart(), one.GetChromEnd(), two.GetChromStart(), two.GetChromEnd())
	expected := 3
	if overlapSize != expected {
		t.Errorf("Error: MatchOverlapLen(%d, %d, %d, %d) %d != %d\n", one.GetChromStart(), one.GetChromEnd(), two.GetChromStart(), two.GetChromEnd(), overlapSize, expected)
	}
}
