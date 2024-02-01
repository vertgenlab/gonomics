package wig

import (
	"github.com/vertgenlab/gonomics/exception"
	"os"
	"testing"

	"github.com/vertgenlab/gonomics/fileio"
)

var SmoothTests = []struct {
	InFile       string
	ChromSizes   string
	OutFile      string
	WindowSize   int
	Missing      float64
	ExpectedFile string
}{
	{InFile: "testdata/unsmooth.wig",
		ChromSizes:   "testdata/smooth.chrom.sizes",
		OutFile:      "testdata/tmp.smooth.wig",
		WindowSize:   5,
		Missing:      -10,
		ExpectedFile: "testdata/expected.smooth.wig"},
}

func TestSmoothMap(t *testing.T) {
	var err error
	var records map[string]Wig
	for _, v := range SmoothTests {
		records = Read(v.InFile, v.ChromSizes, v.Missing)
		records = SmoothMap(records, v.WindowSize, v.Missing)
		Write(v.OutFile, records)
		if !fileio.AreEqual(v.OutFile, v.ExpectedFile) {
			t.Errorf("Error in SmoothSlice. Output not as expected.")
		} else {
			err = os.Remove(v.OutFile)
			exception.PanicOnErr(err)
		}
	}
}
