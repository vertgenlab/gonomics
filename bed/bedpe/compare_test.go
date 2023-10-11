package bedpe

import (
	"testing"

	"github.com/vertgenlab/gonomics/bed"
)

var Equal1 BedPe = BedPe{A: bed.Bed{Chrom: "chr9", ChromStart: 60, ChromEnd: 75}, B: bed.Bed{Chrom: "chr9", ChromStart: 7500, ChromEnd: 7900}}
var Equal2 BedPe = BedPe{A: bed.Bed{Chrom: "chr12", ChromStart: 86500, ChromEnd: 86540}, B: bed.Bed{Chrom: "chr12", ChromStart: 975340, ChromEnd: 975550}}

var EqualTests = []struct {
	a        BedPe
	b        BedPe
	Expected bool
}{
	{Equal1, Equal2, false},
	{Equal2, Equal2, true},
	{Equal1, Equal1, true},
}

func TestEqual(t *testing.T) {
	var actual bool
	for _, v := range EqualTests {
		actual = Equal(v.a, v.b)
		if actual != v.Expected {
			t.Errorf("ERror in bedPe package Equal function.")
		}
	}
}

var AllEqualTests = []struct {
	a        []BedPe
	b        []BedPe
	Expected bool
}{
	{a: []BedPe{Equal1, Equal2, Equal2}, b: []BedPe{Equal1}, Expected: false},
	{a: []BedPe{Equal2, Equal1}, b: []BedPe{Equal2, Equal1}, Expected: true},
	{a: []BedPe{Equal2, Equal2, Equal1, Equal1}, b: []BedPe{Equal2, Equal1, Equal1, Equal1}, Expected: false},
}

func TestAllAreEqual(t *testing.T) {
	var actual bool
	for _, v := range AllEqualTests {
		actual = AllAreEqual(v.a, v.b)
		if actual != v.Expected {
			t.Errorf("Error in bedPe package AllAreEqual function.")
		}
	}
}

var s = []struct {
	truth              []BedPe
	test               []bed.Bed
	expectedFreq       float64
	expectedBedMatches []bed.Bed
}{
	{truth: Read("testdata/statsIn.bedpe"),
		test:               bed.Read("testdata/bedTestIn.bed"),
		expectedFreq:       1.0,
		expectedBedMatches: bed.Read("testdata/expectedMatches.bed")},
}

func TestGeneAssignmentCheck(t *testing.T) {
	var freq float64
	var matches []bed.Bed
	for v := range s {
		freq, matches, _ = GeneAssignmentCheck(s[v].truth, s[v].test)
		bed.Write("testdata/tmp.test.bed", matches)
		if freq != s[v].expectedFreq {
			t.Errorf("frequency found by function did not match expected value. Expected: %f, calculated: %f", s[v].expectedFreq, freq)
		} else if !bed.AllAreEqual(matches, s[v].expectedBedMatches) {
			bed.Write("testdata/tmp.test.bed", matches)
			t.Errorf("Matches found did not match expected bed matches. See testdata/tmp.test.bed.")
		}
	}
}
