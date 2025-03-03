package main

import (
	"os"
	"testing"

	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
)

var LocateCGTest = []struct {
	inFa        string
	chromName   string
	outputTxt   string
	expectedTxt string
}{
	{"testdata/locateCG_test.fa", "chr8", "chr8CGs.txt", "testdata/expected.txt"},
}

func TestLocateCG(t *testing.T) {
	var err error
	var s Settings
	for _, v := range LocateCGTest {
		s = Settings{
			InFa:      v.inFa,
			ChromName: v.chromName,
		}
		locateCG(s)
		if !fileio.AreEqual(v.expectedTxt, v.outputTxt) {
			t.Errorf("Error in locateCG.go, expected.txt != output.txt")
		} else {
			err = os.Remove(v.outputTxt)
			exception.PanicOnErr(err)
		}
	}
}
