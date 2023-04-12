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
}{
	{"testdata/input.fa", "testdata/output.fa", "testdata/expected.fa"},
}

func TestFaInfo(t *testing.T) {
	var err error
	for _, v := range FaInfoTests {
		faInfo(v.inputFile, v.outputFile)
		if !fileio.AreEqual(v.outputFile, v.expectedFile) {
			t.Errorf("Error in faInfo, output and expected do not match.")
		} else {
			err = os.Remove(v.outputFile)
			exception.PanicOnErr(err)
		}
	}
}
