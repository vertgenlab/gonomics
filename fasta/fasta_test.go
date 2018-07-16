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
	data     []Fasta
}{
	{"testdata/testOne.fa", []Fasta{{"apple", seqOneA}, {"banana", seqOneB}, {"carrot", seqOneC}}},
}

func TestRead(t *testing.T) {
	for _, test := range readWriteTests {
		actual, err := Read(test.filename)
		if err != nil {
			t.Errorf("Reading %s gave an error..", test.filename)
		}
		if !AllAreEqual(test.data, actual) {
			t.Errorf("The %s file was not read correctly.", test.filename)
		}
	}
}

func TestWriteAndRead(t *testing.T) {
	var actual []Fasta
	for _, test := range readWriteTests {
		tempFile := test.filename + ".tmp"
		err := Write(tempFile, test.data)
		if err != nil {
			t.Errorf("Error writing %s as a temp fasta file", tempFile)
		}
		actual, err = Read(tempFile)
		if err != nil {
			t.Errorf("Reading %s gave an error", test.filename)
		}
		if !AllAreEqual(test.data, actual) {
			t.Errorf("The %s file was not read correctly.", test.filename)
		}
		err = os.Remove(tempFile)
		if err != nil {
			t.Errorf("Deleting temp file %s gave an error.", tempFile)
		}
	}
}
