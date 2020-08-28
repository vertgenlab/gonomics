package simulate

import (
	"fmt"
	"github.com/vertgenlab/gonomics/dna"
	"testing"
	//"fmt"
)

var GCcontent = 0.42

var RandGeneTests = []struct {
	name   string
	length int
	GC     float64
}{
	{"testingRangGene", 30, GCcontent},
}

func TestRandGene(t *testing.T) {
	for _, test := range RandGeneTests {
		a := RandGene(test.name, test.length, test.GC)
		if len(a[0].Seq) != test.length {
			t.Errorf("expected RandGene to give %v, gave %v", test.length, len(a[0].Seq))
		}
		//fmt.Print(dna.BasesToString(a[0].Seq), "\n")
	}
}

var seq = dna.StringToBases("ATGAGGTCACGATATTAG")
var MutateSeqTests = []struct {
	sequence     []dna.Base
	branchLength float64
}{
	{seq, 1.0}, //branch length of 1 gives higher chance of returning a new base so you can see a difference even with a short sequence
}

func TestMutateSeq(t *testing.T) {
	for _, test := range MutateSeqTests {
		a := MutateSeq(test.sequence, test.branchLength)
		if len(seq) != len(a) {
			t.Errorf("Expected same length sequences. Original: %s \n Ending: %s", seq, dna.BasesToString(a))
		}
		fmt.Print(dna.BasesToString(a), "\n")
	}
}

//tests for functions that MutateSeq is dependent on
//var base = dna.G
//var changeBaseTests = []struct {
//	originalBase dna.Base
//}{
//	{base},
//}
//
//func TestChangeBase(t *testing.T) {
//	for _, test := range changeBaseTests {
//		a := changeBase(test.originalBase)
//		if a == 2 {
//			t.Errorf("Function should have changed base, and didn't.")
//		}
//		fmt.Print(a, "\n")
//	}
//}
//
//var originalBase = dna.G
//var mutateBaseTests = []struct {
//	base         dna.Base
//	branchLength float64
//}{
//	{base, 1.0},
//}

//func TestMutateBase(t *testing.T) {
//	for _, test := range mutateBaseTests {
//		a := mutateBase(test.base, test.branchLength)
//		fmt.Print(a, "\n")
//	}
//}
