package main

import (
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"os"
	"testing"
)

var MergeMultiFaTests = []struct {
	InAFile      string
	InBFile      string
	OutFile      string
	ExpectedFile string
}{
	{InAFile: "testdata/testA.fa",
		InBFile:      "testdata/testB.fa",
		OutFile:      "testdata/tmp.out.fa",
		ExpectedFile: "testdata/expected.out.fa"},
}

func TestMergeMultiFa(t *testing.T) {
	var err error
	for _, v := range MergeMultiFaTests {
		mergeMultiFa(v.InAFile, v.InBFile, v.OutFile)
		if !fileio.AreEqual(v.OutFile, v.ExpectedFile) {
			t.Errorf("Error in mergeMultiFa. Output was not as expected.")
		} else {
			err = os.Remove(v.OutFile)
			exception.PanicOnErr(err)
		}
	}
}
