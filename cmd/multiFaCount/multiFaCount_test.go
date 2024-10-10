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
	either          bool
	secondQueryName string
	thirdQueryName  string
	expectedFile    string
}{
	{"testdata/testInput.fa", "testdata/out.bed", "gibbon", false, false, "", "", "testdata/expected.txt"},
	{"testdata/testInput.fa", "testdata/out2.bed", "orangutan", false, false, "", "", "testdata/expected2.txt"},
	{"testdata/testInput.fa", "testdata/out3.bed", "gibbon", true, false, "orangutan", "", "testdata/expected3.txt"},
	{"testdata/testInput2.fa", "testdata/out4.bed", "gorilla", false, true, "orangutan", "gibbon", "testdata/expected4.txt"},
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
			Either:          v.either,
			SecondQueryName: v.secondQueryName,
			ThirdQueryName:  v.thirdQueryName,
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
