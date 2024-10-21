package main

import (
	"github.com/vertgenlab/gonomics/fileio"
	"os"
	"testing"

	"github.com/vertgenlab/gonomics/exception"
)

var MfaScanTests = []struct {
	inputFile    string
	outputFile   string
	queryName    string
	chrom        string
	expectedFile string
}{
	{"testdata/testInput.fa", "testdata/out.bed", "Human_Chimp_Ancestor", "chr1", "testdata/expected.bed"},
}

func TestMfaScan(t *testing.T) {
	var s Settings
	var err error

	for _, v := range MfaScanTests {

		s = Settings{
			InFile:    v.inputFile,
			OutFile:   v.outputFile,
			QueryName: v.queryName,
			Chrom:     v.chrom,
		}

		multiFaScan(s)

		if !fileio.AreEqual(v.outputFile, v.expectedFile) {
			t.Errorf("Error: output was not as expected.")
		} else {
			err = os.Remove(v.outputFile)
			exception.PanicOnErr(err)
		}
	}
}
