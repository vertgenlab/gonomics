package fasta

import (
	"os"
	"testing"

	"github.com/vertgenlab/gonomics/dna"
)

var seqOneA = dna.StringToBases("ACGTacgTCATCATCATTACTACTAC")
var seqOneB = dna.StringToBases("acgtACGTACGT")
var seqOneC = dna.StringToBases("ACGTACGTACGTT")
var readWriteTests = []struct {
	filename string // input
	data     []Fasta
}{
	{"testdata/testOne.fa", []Fasta{{"apple", seqOneA}, {"banana", seqOneB}, {"carrot", seqOneC}}},
}

func TestRead(t *testing.T) {
	for _, test := range readWriteTests {
		actual := Read(test.filename)
		if !AllAreEqual(test.data, actual) {
			t.Errorf("The %s file was not read correctly.", test.filename)
		}
	}
}

func TestReadToString(t *testing.T) {
	for _, test := range readWriteTests {
		actual := ReadToString(test.filename)
		if len(actual) != len(test.data) {
			t.Error("problem with ReadToString")
		}
		for _, expected := range test.data {
			if dna.BasesToString(expected.Seq) != actual[expected.Name] {
				t.Errorf("problem with ReadToString\n"+
					"expected: %s\n"+
					"actual  : %s\n", dna.BasesToString(expected.Seq), actual[expected.Name])
			}
		}
	}
}

func TestWriteAndRead(t *testing.T) {
	var actual []Fasta
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

func TestToMap(t *testing.T) {
	for _, test := range readWriteTests {
		actual := Read(test.filename)
		faMap := ToMap(actual)
		for key, val := range faMap {
			switch key {
			case test.data[0].Name:
				if dna.CompareSeqsCaseSensitive(val, test.data[0].Seq) != 0 {
					t.Errorf("Problem with ToMap function: seqs do not match\nActual: %s\nExpect: %s\n",
						dna.BasesToString(val), dna.BasesToString(test.data[0].Seq))
				}

			case test.data[1].Name:
				if dna.CompareSeqsCaseSensitive(val, test.data[1].Seq) != 0 {
					t.Errorf("Problem with ToMap function: seqs do not match\nActual: %s\nExpect: %s\n",
						dna.BasesToString(val), dna.BasesToString(test.data[1].Seq))
				}

			case test.data[2].Name:
				if dna.CompareSeqsCaseSensitive(val, test.data[2].Seq) != 0 {
					t.Errorf("Problem with ToMap function: seqs do not match\nActual: %s\nExpect: %s\n",
						dna.BasesToString(val), dna.BasesToString(test.data[2].Seq))
				}

			default:
				t.Errorf("Problem with ToMap: map was larger than expected")
			}
		}
	}
}

var extractTests = []struct {
	fa       Fasta
	start    int
	end      int
	name     string
	expected Fasta
}{
	{fa: Fasta{Name: "chr1", Seq: dna.StringToBases("ATCA")}, start: 0, end: 1, name: "chr1_subset", expected: Fasta{Name: "chr1_subset", Seq: dna.StringToBases("A")}},
	{fa: Fasta{Name: "chr1", Seq: dna.StringToBases("ATCA")}, start: 3, end: 4, name: "chr1_lastBase", expected: Fasta{Name: "chr1_lastBase", Seq: dna.StringToBases("A")}}, // last base
}

func TestExtract(t *testing.T) {
	var actual Fasta
	for _, v := range extractTests {
		actual = Extract(v.fa, v.start, v.end, v.name)
		if !IsEqual(actual, v.expected) {
			t.Errorf("Error in testExtract. Actual: %v. Expected: %v.\n", actual, v.expected)
		}
	}
}
