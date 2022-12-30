package motif

import (
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"os"
	"testing"
)

var ReadWritePfmTests = []struct {
	InFile       string
	OutFile      string
	ExpectedFile string
	Type string
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
		records = Read(v.InFile, v.Type)
		Write(v.OutFile, records)
		if !fileio.AreEqual(v.OutFile, v.ExpectedFile) {
			t.Errorf("Error in Read/Write PositionMatrix. Output not as expected.")
		} else {
			err = os.Remove(v.OutFile)
			exception.PanicOnErr(err)
		}
	}
}
