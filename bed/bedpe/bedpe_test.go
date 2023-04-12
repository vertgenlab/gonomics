package bedpe

import (
	"os"
	"testing"

	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/exception"
)

var b1 BedPe = BedPe{
	A: bed.Bed{
		Chrom:             "chr2",
		ChromStart:        200,
		ChromEnd:          240,
		Name:              "First",
		Score:             40,
		Strand:            '+',
		FieldsInitialized: 9,
	},
	B: bed.Bed{
		Chrom:             "chr2",
		ChromStart:        4500,
		ChromEnd:          4600,
		Name:              "First",
		Score:             40,
		Strand:            '-',
		FieldsInitialized: 9,
	},
}

var b2 BedPe = BedPe{
	A: bed.Bed{
		Chrom:             "chr9",
		ChromStart:        3400,
		ChromEnd:          3500,
		Name:              "Second",
		Score:             2,
		Strand:            '-',
		FieldsInitialized: 9,
	},
	B: bed.Bed{
		Chrom:             "chr9",
		ChromStart:        6700,
		ChromEnd:          6780,
		Name:              "Second",
		Score:             2,
		Strand:            '-',
		FieldsInitialized: 9,
	},
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
