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
	{"testingRandGene", 102, GCcontent},
}

func TestRandGene(t *testing.T) {
	for _, test := range RandGeneTests {
		a := RandGene(test.name, test.length, test.GC)
		if len(a[0].Seq) != test.length {
			t.Errorf("expected RandGene to give %v, gave %v", test.length, len(a[0].Seq))
		}
		//fmt.Print(dna.BasesToString(a[0].Seq), "\n")
		//fasta.Write("outputRandGene.fasta", a)
	}
}

var MutateSeqTests = []struct {
	sequence     string
	branchLength float64
	gp           string
}{
	{"testdata/longDebug.fasta", 0.5, "testdata/longDebug.gp"}, //branch length of 1 gives higher chance of returning a new base so you can see a difference even with a short sequence
}

//func TestBaseConversions(t *testing.T) {
//	for _, test := range MutateSeqTests {
//		seq := fasta.Read(test.sequence)
//		bases := seq[0].Seq
//		intermediate := BasesToBaseExt(bases)
//		answer := BaseExtToBases(intermediate)
//		for i := 0; i < len(intermediate); i++ {
//			if intermediate[i].Base != answer[i] {
//				log.Fatal("bases are not in the right order")
//			}
//		}
//		if len(answer) != len(bases) {
//			log.Fatal("Test Base Conversion: answer not written properly")
//		}
//	}
//}

//func TestCodonExtConversions(t *testing.T) {
//	for _, test := range MutateSeqTests {
//		fasta := fasta.Read(test.sequence)
//		seq := BasesToBaseExt(fasta[0].Seq)
//		gene := genePred.Read(test.gp)
//		a := CodonExtsToCodons(CreateCodons(seq, gene[0], 0))
//		log.Print(a)
//		b := CodonExtToBaseExt(CreateCodons(seq, gene[0], 0))
//		log.Print(b)
//	}
//}

func TestMutateGene(t *testing.T) {
	for _, test := range MutateSeqTests {
		seq := fasta.Read(test.sequence)
		bases := seq[0].Seq
		a := MutateGene(bases, test.branchLength, test.gp, true)
		fmt.Printf("a: %s\n", dna.BasesToString(a))
		if len(bases) != len(a) {
			t.Errorf("Expected same length sequences. Original: %v \n Ending: %v", len(bases), len(a))
		}
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
//	position	int
//}{
//	{originalBase, 1.0, 0},
//}
//
//func TestMutateBase(t *testing.T) {
//	for _, test := range mutateBaseTests {
//		a := mutateBase(test.base, test.branchLength, test.position)
//		fmt.Print(a, "\n")
//	}
//}
