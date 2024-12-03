package lift

import (
	"github.com/vertgenlab/gonomics/axt"
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

func TestLiftCoordinatesWithAxt(t *testing.T) {
	var outBed bed.Bed = bed.Bed{FieldsInitialized: 3}
	axtRecords := axt.Read("testdata/in.axt")
	bedsToLift := []bed.Bed{{Chrom: "chr1", ChromStart: 100, ChromEnd: 110}, {Chrom: "chr1", ChromStart: 105, ChromEnd: 110}}
	expectedBed := []bed.Bed{{Chrom: "chr2", ChromStart: 200, ChromEnd: 210, FieldsInitialized: 3}, {Chrom: "chr2", ChromStart: 789, ChromEnd: 796, FieldsInitialized: 3}}
	for i := range bedsToLift {
		outBed.Chrom, outBed.ChromStart, outBed.ChromEnd = LiftCoordinatesWithAxt(axtRecords[i], bedsToLift[i], 1000)
		if !bed.Equal(expectedBed[i], outBed) {
			t.Errorf("Error in LiftCoordinatesWithAxt. Expected: %v, output: %v\n", expectedBed[i], outBed)
		}
	}
}
