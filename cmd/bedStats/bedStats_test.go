package main

import (
	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/bed/bedpe"
	"os"
	"testing"
)

var s = []struct {
	truth              []bedpe.BedPe
	test               []bed.Bed
	expectedFreq       float64
	expectedBedMatches []bed.Bed
	expectedNonMatches []bed.Bed
}{
	{truth: bedpe.Read("testdata/statsIn.bedpe"),
		test:               bed.Read("testdata/bedTestIn.bed"),
		expectedFreq:       1.0,
		expectedBedMatches: bed.Read("testdata/expectedMatches.bed"),
		expectedNonMatches: bed.Read("testdata/expectedNonMatches.bed")},
}

func TestGeneAssignmentCheck(t *testing.T) {
	var freq float64
	var matches, nonMatches []bed.Bed
	for v := range s {
		freq, matches, nonMatches = GeneAssignmentCheck(s[v].truth, s[v].test)
		bed.Write("testdata/tmp.test.bed", matches)
		if freq != s[v].expectedFreq {
			//t.Errorf("frequency found by function did not match expected value. Expected: %f, calculated: %f", s[v].expectedFreq, freq)
		} else if !bed.AllAreEqual(matches, s[v].expectedBedMatches) {
			bed.Write("testdata/tmp.test.bed", matches)
			t.Errorf("Matches found did not match expected bed matches. See testdata/tmp.test.bed.")
		} else if !bed.AllAreEqual(nonMatches, s[v].expectedNonMatches) {
			bed.Write("testdata/tmp.test.bed", nonMatches)
			t.Errorf("Unexpected values for nonMatching beds. See testdata/tmp.test.bed.")
		} else {
			os.Remove("testdata/tmp.test.bed")
		}
	}
}
