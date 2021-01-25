package fasta

import (
	"github.com/vertgenlab/gonomics/dna"
	"os"
	"testing"
)

var seqOneA, _ = dna.StringToBases("ACGTacgTCATCATCATTACTACTAC")
var seqOneB, _ = dna.StringToBases("acgtACGTACGT")
var seqOneC, _ = dna.StringToBases("ACGTACGTACGTT")
var readWriteTests = []struct {
	filename string // input
	data     []*Fasta
}{
	{"testdata/testOne.fa", []*Fasta{{"apple", seqOneA}, {"banana", seqOneB}, {"carrot", seqOneC}}},
}

func TestRead(t *testing.T) {
	for _, test := range readWriteTests {
		actual := Read(test.filename)
		if !AllAreEqual(test.data, actual) {
			t.Errorf("The %s file was not read correctly.", test.filename)
		}
	}
}

func TestWriteAndRead(t *testing.T) {
	var actual []*Fasta
	for _, test := range readWriteTests {
		tempFile := test.filename + ".tmp"
		Write(tempFile, test.data)
		actual = Read(tempFile)
		if !AllAreEqual(test.data, actual) {
			t.Errorf("The %s file was not read correctly.", test.filename)
		}
		err := os.Remove(tempFile)
		if err != nil {
			t.Errorf("Deleting temp file %s gave an error.", tempFile)
		}
	}
}
