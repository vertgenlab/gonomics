package main

import (
	"github.com/vertgenlab/gonomics/fileio"
	"os"
	"testing"

	"github.com/vertgenlab/gonomics/exception"
)

var MfaCountTests = []struct {
	inputFile       string
	outputFile      string
	queryName       string
	both            bool
	secondQueryName string
	expectedFile    string
}{
	{"testdata/testInput.fa", "testdata/out.bed", "gibbon", false, "", "testdata/expected.txt"},
	{"testdata/testInput.fa", "testdata/out2.bed", "orangutan", false, "", "testdata/expected2.txt"},
	{"testdata/testInput.fa", "testdata/out3.bed", "gibbon", true, "orangutan", "testdata/expected3.txt"},
}

func TestMfaCount(t *testing.T) {
	var s Settings
	var err error

	for _, v := range MfaCountTests {

		s = Settings{
			InFile:          v.inputFile,
			OutFile:         v.outputFile,
			QueryName:       v.queryName,
			Both:            v.both,
			SecondQueryName: v.secondQueryName,
		}

		multiFaCount(s)

		if !fileio.AreEqual(v.outputFile, v.expectedFile) {
			t.Errorf("Error: output was not as expected.")
		} else {
			err = os.Remove(v.outputFile)
			exception.PanicOnErr(err)
		}
	}
}
