package main

import (
	"os"
	"testing"

	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
)

var FormatIdeogramTests = []struct {
	InFile       string
	OutFile      string
	ExpectedFile string
	NoScore      bool
}{
	{"testdata/test.bed", "testdata/tmp.txt", "testdata/expected.Score.txt", false},
	{"testdata/test.bed", "testdata/tmp.txt", "testdata/expected.NoScore.txt", true},
}

func TestFormatIdeogram(t *testing.T) {
	var err error
	for _, v := range FormatIdeogramTests {
		formatIdeogram(v.InFile, v.OutFile, v.NoScore)
		if !fileio.AreEqual(v.OutFile, v.ExpectedFile) {
			t.Errorf("Error in formatIdeogram. Output does not match expected.")
		} else {
			err = os.Remove(v.OutFile)
			exception.PanicOnErr(err)
		}
	}
}
