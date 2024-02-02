package main

import (
	"os"
	"testing"

	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
)

var FaInfoTests = []struct {
	inputFile    string
	outputFile   string
	expectedFile string
	gcContent    bool
}{
	{"testdata/input.fa", "testdata/output.fa", "testdata/expected.fa", true},
}

func TestFaInfo(t *testing.T) {
	var err error
	for _, v := range FaInfoTests {
		faInfo(v.inputFile, v.outputFile, v.gcContent)
		if !fileio.AreEqual(v.outputFile, v.expectedFile) {
			t.Errorf("Error in faInfo, output and expected do not match.")
		} else {
			err = os.Remove(v.outputFile)
			exception.PanicOnErr(err)
		}
	}
}

//TODO: update testing files for options
