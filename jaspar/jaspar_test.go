package jaspar

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
}{
	{"testdata/jaspar.vertebrate.txt",
		"testdata/tmp.jaspar.txt",
		"testdata/expected.jaspar.txt"},
}

func TestReadAndWritePfm(t *testing.T) {
	var err error
	var records []Pfm
	for _, v := range ReadWritePfmTests {
		records = ReadPfm(v.InFile)
		WritePfmSlice(v.OutFile, records)
		if !fileio.AreEqual(v.OutFile, v.ExpectedFile) {
			t.Errorf("Error in ReadWritePfm. Output not as expected.")
		} else {
			err = os.Remove(v.OutFile)
			exception.PanicOnErr(err)
		}
	}
}
