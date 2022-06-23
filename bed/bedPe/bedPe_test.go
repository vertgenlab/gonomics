package bedpe

import (
	"github.com/vertgenlab/gonomics/exception"
	"os"
	"testing"
)

var b1 BedPe = BedPe{ChromA: "chr2",
	ChromStartA:       200,
	ChromEndA:         240,
	ChromB:            "chr2",
	ChromStartB:       4500,
	ChromEndB:         4600,
	Name:              "First",
	Score:             40,
	StrandA:           '+',
	StrandB:           '-',
	FieldsInitialized: 9,
}

var b2 BedPe = BedPe{ChromA: "chr9",
	ChromStartA: 3400,
	ChromEndA: 3500,
	ChromB: "chr9",
	ChromStartB: 6700,
	ChromEndB: 6780,
	Name: "Second",
	Score: 2,
	StrandA: '-',
	StrandB: '-',
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
