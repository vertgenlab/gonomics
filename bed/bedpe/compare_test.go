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
