package bedPe

import (
	"testing"
)

var Equal1 BedPe = BedPe{Chrom1: "chr9", ChromStart1: 60, ChromEnd1: 75, Chrom2: "chr9", ChromStart2: 7500, ChromEnd2: 7900}
var Equal2 BedPe = BedPe{Chrom1: "chr12", ChromStart1: 86500, ChromEnd1: 86540, Chrom2: "chr12", ChromStart2: 975340, ChromEnd2: 975550}

var EqualTests = []struct {
	a BedPe
	b BedPe
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
	a []BedPe
	b []BedPe
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
