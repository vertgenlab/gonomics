package main

import (
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"os"
	"testing"
)

var ScCountTests = []struct {
	InFile string
	OutFile string
	ExpectedFile string
	GeneFile string
}{
	{"testdata/test.sam", "testdata/out.tsv", "testdata/expected.tsv", "testdata/test.gtf"},
}

func TestScCount(t *testing.T) {
	var err error
	var s Settings
	for _, v := range ScCountTests {
		s = Settings{
			InFile: v.InFile,
			OutFile: v.OutFile,
			GeneFile: v.GeneFile,
		}
		scCount(s)
		if !fileio.AreEqual(v.OutFile, v.ExpectedFile) {
			t.Errorf("Error in scCount.")
		} else {
			err = os.Remove(v.OutFile)
			exception.PanicOnErr(err)
		}
	}
}