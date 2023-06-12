package motif

import (
	"os"
	"testing"

	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
)

var ReadWritePfmTests = []struct {
	InFile       string
	OutFile      string
	ExpectedFile string
	Type         string
}{
	{"testdata/jaspar.vertebrate.txt",
		"testdata/tmp.jaspar.txt",
		"testdata/expected.jaspar.txt",
		"Frequency"},
}

func TestReadAndWrite(t *testing.T) {
	var err error
	var records []PositionMatrix
	for _, v := range ReadWritePfmTests {
		records = ReadJaspar(v.InFile, v.Type)
		WriteJaspar(v.OutFile, records)
		if !fileio.AreEqual(v.OutFile, v.ExpectedFile) {
			t.Errorf("Error in Read/Write PositionMatrix. Output not as expected.")
		} else {
			err = os.Remove(v.OutFile)
			exception.PanicOnErr(err)
		}
	}
}
