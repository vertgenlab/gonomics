package simulate

import (
	"fmt"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/fasta"
	"testing"
	//"fmt"
)

var GCcontent = 0.42

var RandGeneTests = []struct {
	name   string
	length int
	GC     float64
}{
	{"testingRandGene", 99, GCcontent},
}

func TestRandGene(t *testing.T) {
	for _, test := range RandGeneTests {
		a := RandGene(test.name, test.length, test.GC)
		if len(a[0].Seq) != test.length {
			t.Errorf("expected RandGene to give %v, gave %v", test.length, len(a[0].Seq))
		}
		fmt.Print(dna.BasesToString(a[0].Seq), "\n")
	}
}

var MutateSeqTests = []struct {
	sequence     string
	branchLength float64
	gtf          string
}{
	{"debug.fasta", 0.5, "debug.gtf"}, //branch length of 1 gives higher chance of returning a new base so you can see a difference even with a short sequence
}

func TestMutateSeq(t *testing.T) {
	for _, test := range MutateSeqTests {
		seq := fasta.Read(test.sequence)
		bases := seq[0].Seq
		a := MutateSeq(bases, test.branchLength, test.gtf)
		if len(bases) != len(a) {
			t.Errorf("Expected same length sequences. Original: %s \n Ending: %s", dna.BasesToString(bases), dna.BasesToString(a))
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
