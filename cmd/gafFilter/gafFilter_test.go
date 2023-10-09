package main

import (
	"os"
	"testing"

	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
)

var GafFilterTests = []struct {
	InFile       string
	OutFile      string
	ExpectedFile string
	RemoveNot    bool
}{
	{InFile: "testdata/test.gaf",
		OutFile:      "testdata/tmp.gaf",
		ExpectedFile: "testdata/expected.gaf",
		RemoveNot:    true,
	},
}

func TestGafFilter(t *testing.T) {
	var err error
	var s Settings
	for _, v := range GafFilterTests {
		s = Settings{
			InFile:    v.InFile,
			OutFile:   v.OutFile,
			RemoveNot: v.RemoveNot,
		}
		gafFilter(s)
		if !fileio.AreEqual(v.OutFile, v.ExpectedFile) {
			t.Errorf("Error in gafFilter.")
		} else {
			err = os.Remove(v.OutFile)
			exception.PanicOnErr(err)
		}
	}
}
