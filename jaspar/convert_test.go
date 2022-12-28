package jaspar

import (
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"os"
	"testing"
)

var PfmToPpmTests = []struct {
	PfmFile      string
	OutputFile   string
	ExpectedFile string
	Pseudocount  float64
}{
	{"testdata/expected.jaspar.txt",
		"testdata/tmp.Ppm.txt",
		"testdata/expected.Ppm.txt",
		0.1},
}

func TestPfmSliceToPpmSlice(t *testing.T) {
	var err error
	var records []Pfm
	var answer []Ppm
	for _, v := range PfmToPpmTests {
		records = ReadPfm(v.PfmFile)
		answer = PfmSliceToPpmSlice(records, v.Pseudocount)
		WritePpmSlice(v.OutputFile, answer)
		if !fileio.AreEqual(v.OutputFile, v.ExpectedFile) {
			t.Errorf("Error in PfmSliceToPpmSlice. Output was not as expected.")
		} else {
			err = os.Remove(v.OutputFile)
			exception.PanicOnErr(err)
		}
	}
}
