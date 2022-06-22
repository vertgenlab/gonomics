package bedPe

import (
	"github.com/vertgenlab/gonomics/exception"
	"os"
	"testing"
)

var b1 BedPe = BedPe{Chrom1: "chr2",
	ChromStart1:       200,
	ChromEnd1:         240,
	Chrom2:            "chr2",
	ChromStart2:       4500,
	ChromEnd2:         4600,
	Name:              "First",
	Score:             40,
	Strand1:           '+',
	Strand2:           '-',
	FieldsInitialized: 9,
}

var b2 BedPe = BedPe{Chrom1: "chr9",
	ChromStart1: 3400,
	ChromEnd1: 3500,
	Chrom2: "chr9",
	ChromStart2: 6700,
	ChromEnd2: 6780,
	Name: "Second",
	Score: 2,
	Strand1: '-',
	Strand2: '-',
	FieldsInitialized: 9,
	}

var beds []BedPe = []BedPe{b1, b2}

var readWriteTests = []struct {
	filename string
	data     []BedPe
}{
	{"testdata/BedPeFileTest.bedpe", beds},
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
	var actual []BedPe
	var err error
	for _, test := range readWriteTests {
		tempFile := test.filename + ".tmp"
		Write(tempFile, test.data)
		actual = Read(test.filename)
		if !AllAreEqual(test.data, actual) {
			t.Errorf("The %s file was not read correctly.", test.filename)
		}
		err = os.Remove(tempFile)
		exception.PanicOnErr(err)
	}

}
