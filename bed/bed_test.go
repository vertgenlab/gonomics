package bed

import (
	"os"
	"testing"
)

var b1 Bed = Bed{Name: "chr1", Start: 100, End: 200}
var b2 Bed = Bed{Name: "chr2", Start: 400, End: 900}
var b3 Bed = Bed{Name: "chr3", Start: 945, End: 1000}
var beds []*Bed = []*Bed{&b1, &b2, &b3}

var readWriteTests = []struct {
	filename string
	data     []*Bed
}{
	{"testdata/bedFileTest.bed", beds},
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
	var actual []*Bed
	for _, test := range readWriteTests {
		tempFile := test.filename + ".tmp"
		err := Write(tempFile, test.data)
		if err != nil {
			t.Errorf("Error writing %s as a temp bed file", tempFile)
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
