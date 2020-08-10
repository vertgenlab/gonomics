package simulate

import (
	"testing"
	"github.com/vertgenlab/gonomics/dna"
	"log"
	"fmt"
)


var RandGeneTests = []struct {
	name string
	length int
}{
	{"testingRangGene", 30},
}

func TestRandGene(t *testing.T) {
	for _, test := range RandGeneTests {
		a := RandGene(test.name, test.length)
		if len(a[0].Seq) != test.length {
			t.Errorf("expected RandGene to give %v, gave %v", test.length, len(a[0].Seq))
		}
	}
}

var seq = dna.StringToBases("ATGAGGTCACGATATTAG")
var MutateSeqTests = []struct {
	sequence []dna.Base
	branchLength float64
}{
	{seq, 1.0},//branch length of 1 gives higher chance of returning a new base so you can see a difference even with a short sequence
}

func TestMutateSeq(t *testing.T) {
	for _, test := range MutateSeqTests {
		a := MutateSeq(test.sequence, test.branchLength)
		if len(seq) != len(a) {
			t.Errorf("Expected same length sequences. Original: %s \n Ending: %s", seq, dna.BasesToString(a))
		}
		log.Printf(dna.BasesToString(a))
	}
}

var base = dna.G
var changeBaseTests = []struct {
	originalBase dna.Base
} {
	{base},
}

func TestChangeBase(t *testing.T) {
	for _, test := range changeBaseTests {
		a := changeBase(test.originalBase)
		if a == 2 {
			t.Errorf("Function should have changed base, and didn't.")
		}
		fmt.Print(a, "\n")
	}
}

var originalBase = dna.G
var mutateBaseTests = []struct {
	base dna.Base
	branchLength float64
} {
	{base, 1.0},
}

func TestMutateBase(t *testing.T) {
	for _, test := range mutateBaseTests {
		a := mutateBase(test.base, test.branchLength)
		fmt.Print(a, "\n")
	}
}