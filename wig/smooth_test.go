package wig

import (
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"os"
	"testing"
)

var SmoothTests = []struct {
	InFile       string
	OutFile      string
	WindowSize   int
	Missing float64
	ExpectedFile string
}{
	{InFile: "testdata/unsmooth.wig",
		OutFile:      "testdata/tmp.smooth.wig",
		WindowSize:   5,
		Missing: -10,
		ExpectedFile: "testdata/expected.smooth.wig"},
}

func TestSmoothSlice(t *testing.T) {
	var err error
	var records []Wig
	for _, v := range SmoothTests {
		records = Read(v.InFile)
		records = SmoothSlice(records, v.WindowSize, v.Missing)
		Write(v.OutFile, records)
		if !fileio.AreEqual(v.OutFile, v.ExpectedFile) {
			t.Errorf("Error in SmoothSlice. Output not as expected.")
		} else {
			err = os.Remove(v.OutFile)
			exception.PanicOnErr(err)
		}
	}
}
