package main

import (
	"os"
	"testing"

	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
)

var DotToNewickTests = []struct {
	InFile       string
	OutFile      string
	ExpectedFile string
	Verbose      bool
}{
	{"testdata/primate.dot", "testdata/tmp.nh", "testdata/expected.nh", false},
}

func TestDotToNewick(t *testing.T) {
	var err error
	for _, v := range DotToNewickTests {
		dotToNewick(v.InFile, v.OutFile, v.Verbose)
		if !fileio.AreEqual(v.OutFile, v.ExpectedFile) {
			t.Errorf("Error in dotToNewick. Output did not match expected.")
		} else {
			err = os.Remove(v.OutFile)
			exception.PanicOnErr(err)
		}
	}
}
