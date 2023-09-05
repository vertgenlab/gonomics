package fasta

import (
	"testing"

	"github.com/vertgenlab/gonomics/dna"
)

var seqTwoA = dna.StringToBases("ACGTacgTCATCATCATTACTACTAC")
var seqTwoB = dna.StringToBases("acgtACGTACGT")
var seqTwoC = dna.StringToBases("ACGTACGTACGTT")
var allEqualTests = []struct {
	dataOne  []Fasta
	dataTwo  []Fasta
	expected bool
}{
	{[]Fasta{{"apple", seqTwoA}, {"banana", seqTwoB}, {"carrot", seqTwoC}}, []Fasta{{"apple", seqTwoA}, {"banana", seqTwoB}, {"carrot", seqTwoC}}, true},
	{[]Fasta{{"apple", seqTwoA}, {"bananaa", seqTwoB}, {"carrot", seqTwoC}}, []Fasta{{"apple", seqTwoA}, {"banana", seqTwoB}, {"carrot", seqTwoC}}, false},
	{[]Fasta{{"apple", seqTwoA}, {"banana", seqTwoB}, {"carrot", seqTwoC}}, []Fasta{{"apple", seqTwoA}, {"banana", seqTwoB}, {"carrot", seqTwoB}}, false},
}

func TestAllEqual(t *testing.T) {
	for _, test := range allEqualTests {
		actual := AllAreEqual(test.dataOne, test.dataTwo)
		if actual != test.expected {
			t.Errorf("Compared %v and %v; expected %t and got %t.", test.dataOne, test.dataTwo, test.expected, actual)
		}
	}
}

var sortNameTests = []struct {
	input    []Fasta
	expected []Fasta
}{
	{[]Fasta{{"apple", seqTwoA}, {"banana", seqTwoB}, {"carrot", seqTwoC}}, []Fasta{{"apple", seqTwoA}, {"banana", seqTwoB}, {"carrot", seqTwoC}}},
	{[]Fasta{{"banana", seqTwoB}, {"apple", seqTwoA}, {"carrot", seqTwoC}}, []Fasta{{"apple", seqTwoA}, {"banana", seqTwoB}, {"carrot", seqTwoC}}},
	{[]Fasta{{"carrot", seqTwoC}, {"apple", seqTwoA}, {"banana", seqTwoB}}, []Fasta{{"apple", seqTwoA}, {"banana", seqTwoB}, {"carrot", seqTwoC}}},
}

func TestSortName(t *testing.T) {
	for _, test := range sortNameTests {
		SortByName(test.input)
		if !AllAreEqual(test.input, test.expected) {
			t.Errorf("Got %v when sorted by name, but expected %v.", test.input, test.expected)
		}
	}
}

var sortSeqTests = []struct {
	input    []Fasta
	expected []Fasta
}{
	{[]Fasta{{"apple", seqTwoA}, {"banana", seqTwoB}, {"carrot", seqTwoC}}, []Fasta{{"banana", seqTwoB}, {"carrot", seqTwoC}, {"apple", seqTwoA}}},
	{[]Fasta{{"banana", seqTwoB}, {"apple", seqTwoA}, {"carrot", seqTwoC}}, []Fasta{{"banana", seqTwoB}, {"carrot", seqTwoC}, {"apple", seqTwoA}}},
	{[]Fasta{{"carrot", seqTwoC}, {"apple", seqTwoA}, {"banana", seqTwoB}}, []Fasta{{"banana", seqTwoB}, {"carrot", seqTwoC}, {"apple", seqTwoA}}},
}

func TestSortSeq(t *testing.T) {
	for _, test := range sortSeqTests {
		SortBySeq(test.input)
		if !AllAreEqual(test.input, test.expected) {
			t.Errorf("Got %v when sorted by sequence, but expected %v.", test.input, test.expected)
		}
	}
}
