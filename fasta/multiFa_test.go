package fasta

import (
	"github.com/vertgenlab/gonomics/dna"
	"testing"
)

var seqThreeA = dna.StringToBases("ACGTacgTCATCATCATTACTACTAC")
var seqThreeB = dna.StringToBases("acgtACGTACGT")
var seqThreeC = dna.StringToBases("ACGTACGTACGTT")
var rcSeqThreeA = dna.StringToBases("GTAGTAGTAATGATGATGAcgtACGT")
var rcSeqThreeB = dna.StringToBases("ACGTACGTacgt")
var rcSeqThreeC = dna.StringToBases("AACGTACGTACGT")
var allRevCompTests = []struct {
	input    []*Fasta
	expected []*Fasta
}{
	{[]*Fasta{{"apple", seqThreeA}, {"banana", seqThreeB}, {"carrot", seqThreeC}}, []*Fasta{{"apple", rcSeqThreeA}, {"banana", rcSeqThreeB}, {"carrot", rcSeqThreeC}}},
}
var segSite1 = dna.StringToBases("AATCCTATTCA")
var segSite2 = dna.StringToBases("AATCCTAATCA")
var segSite3 = dna.StringToBases("AATCCTATTCG")
var expectedSegSite1 = dna.StringToBases("TA")
var expectedSegSite2 = dna.StringToBases("AA")
var expectedSegSite3 = dna.StringToBases("TG")
var segSiteInput = []*Fasta{{"apple", segSite1}, {"banana", segSite2}, {"carrot", segSite3}}
var segSiteExpected = []*Fasta{{"apple", expectedSegSite1}, {"banana", expectedSegSite2}, {"carrot", expectedSegSite3}}

func TestReverseComplement(t *testing.T) {
	for _, test := range allRevCompTests {
		ReverseComplementAll(test.input)
		if !AllAreEqual(test.input, test.expected) {
			t.Errorf("Expected reverse complement to give %v, but got %v.", test.input, test.expected)
		}
	}
}

func TestSegregatingSites(t *testing.T) {
	input := SegregatingSites(segSiteInput)
	if !AllAreEqual(input, segSiteExpected) {
		t.Errorf("Expected segregating sites to give %v columns, but got %v.", len(segSiteExpected[0].Seq), len(input[0].Seq))
	}
}

/*
func TestChangeName(t *testing.T) {
	fa := Read("testdata/testOne.fa")
	newFa := Read("testdata/testOne.fa")
	ChangePrefix(newFa, "Library")
	for i := 0; i < len(newFa); i++ {
		fmt.Printf("Name of fasta record: %s\n", newFa[i].Name)
		if newFa[i].Name == fa[i].Name {
			t.Errorf("Fasta record change was not successful")
		}
	}
}*/
