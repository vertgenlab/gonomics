package fasta

import (
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"os"
	"testing"
)

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

var ReverseTests = []struct {
	input string
	output string
	expected string
}{
	{"testdata/testOne.fa", "testdata/rev.testOne.fa", "testdata/expected.rev.testOne.fa"},
}

func TestReverse(t *testing.T) {
	var records []Fasta
	var err error
	for _, v := range ReverseTests {
		records = Read(v.input)
		for i := range records {
			Reverse(records[i])
		}
		Write(v.output, records)
		if !fileio.AreEqual(v.output, v.expected) {
			t.Errorf("Error in fasta/modify.go. Output of Reverse function was not as expected.")
		} else {
			err = os.Remove(v.output)
			exception.PanicOnErr(err)
		}
	}
}
