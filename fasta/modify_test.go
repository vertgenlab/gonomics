package fasta

import (
	"github.com/vertgenlab/gonomics/dna"
	"testing"
)

var seqThreeA, _ = dna.StringToBases("ACGTacgTCATCATCATTACTACTAC")
var seqThreeB, _ = dna.StringToBases("acgtACGTACGT")
var seqThreeC, _ = dna.StringToBases("ACGTACGTACGTT")
var rcSeqThreeA, _ = dna.StringToBases("GTAGTAGTAATGATGATGAcgtACGT")
var rcSeqThreeB, _ = dna.StringToBases("ACGTACGTacgt")
var rcSeqThreeC, _ = dna.StringToBases("AACGTACGTACGT")
var allRevCompTests = []struct {
	input    []Fasta
	expected []Fasta
}{
	{[]Fasta{{"apple", seqThreeA}, {"banana", seqThreeB}, {"carrot", seqThreeC}}, []Fasta{{"apple", rcSeqThreeA}, {"banana", rcSeqThreeB}, {"carrot", rcSeqThreeC}}},
}

func TestReverseComplement(t *testing.T) {
	for _, test := range allRevCompTests {
		ReverseComplementAll(test.input)
		if !AllAreEqual(test.input, test.expected) {
			t.Errorf("Expected reverse complement to give %v, but got %v.", test.input, test.expected)
		}
	}
}
