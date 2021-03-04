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
	input    []Fasta
	expected []Fasta
}{
	{[]Fasta{{"apple", seqThreeA}, {"banana", seqThreeB}, {"carrot", seqThreeC}}, []Fasta{{"apple", rcSeqThreeA}, {"banana", rcSeqThreeB}, {"carrot", rcSeqThreeC}}},
}

func TestRemove(t *testing.T) {
	for _, test := range readWriteTests {
		rmCopy := test.data[1]
		expected := []Fasta{test.data[0], test.data[2]}
		actual := Remove(test.data, 1)
		if !AllAreEqual(expected, actual) {
			t.Errorf("The %s file was not read correctly.", test.filename)
		}
		// add back so other tests are not broken
		test.data[1] = rmCopy
	}
}

func TestReverseComplement(t *testing.T) {
	for _, test := range allRevCompTests {
		ReverseComplementAll(test.input)
		if !AllAreEqual(test.input, test.expected) {
			t.Errorf("Expected reverse complement to give %v, but got %v.", test.input, test.expected)
		}
	}
}
