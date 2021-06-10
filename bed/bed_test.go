package bed

import (
	"os"
	"testing"
)

var b1 Bed = Bed{Chrom: "chr1", ChromStart: 100, ChromEnd: 200, Name: "First", Score: 1, Strand: '+'}
var b2 Bed = Bed{Chrom: "chr2", ChromStart: 400, ChromEnd: 900, Name: "Second", Score: 5, Strand: '-'}
var b3 Bed = Bed{Chrom: "chr3", ChromStart: 945, ChromEnd: 1000, Name: "Third", Score: 10, Strand: '.'}
var beds []*Bed = []*Bed{&b1, &b2, &b3}

var readWriteTests = []struct {
	filename string
	data     []*Bed
}{
	{"testdata/bedFileTest.bed", beds},
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
	var actual []*Bed
	for _, test := range readWriteTests {
		tempFile := test.filename + ".tmp"
		Write(tempFile, test.data, 3)
		actual = Read(tempFile)
		if !AllAreEqual(test.data, actual) {
			t.Errorf("The %s file was not written and read correctly.", test.filename)
		}
		err := os.Remove(tempFile)
		if err != nil {
			t.Errorf("Deleting temp file %s gave an error.", tempFile)
		}
	}
}
