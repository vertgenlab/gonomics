package maf

import (
	"github.com/vertgenlab/gonomics/fasta"
	"testing"
)

var toFastaTests = []struct {
	inputMaf     string
	inputFa      string
	inputSpecies []string
	expected     string
}{
	{"testdata/toFastaTest.maf", "testdata/toFastaTest.fa", []string{"hg38", "panPan2", "panTro6", "gorGor5", "ponAbe3"}, "testdata/toFastaTest.mfa"},
}

func TestToFasta(t *testing.T) {
	for _, test := range toFastaTests {
		inMaf := Read(test.inputMaf)
		inFa := fasta.Read(test.inputFa)
		actual := ToFasta(inMaf, inFa[0], test.inputSpecies)
		expected := fasta.Read(test.expected)
		if !fasta.AllAreEqual(actual, expected) {
			t.Errorf("Output is not as expected by maf.ToFasta")
		}
	}
}
