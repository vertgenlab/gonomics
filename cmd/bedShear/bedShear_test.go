package main

import (
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"os"
	"testing"
)

var BedShearTests = []struct {
	InFile       string
	OutFile      string
	ExpectedFile string
	FragmentSize int
}{
	{InFile: "testdata/test.bed",
		OutFile:      "testdata/tmp.1.bed",
		ExpectedFile: "testdata/expected.1.bed",
		FragmentSize: 1,
	},
	{InFile: "testdata/test.bed",
		OutFile:      "testdata/tmp.7.bed",
		ExpectedFile: "testdata/expected.7.bed",
		FragmentSize: 7,
	},
}

func TestBedShear(t *testing.T) {
	var err error
	for _, v := range BedShearTests {
		bedShear(v.InFile, v.OutFile, v.FragmentSize)
		if !fileio.AreEqual(v.OutFile, v.ExpectedFile) {
			t.Errorf("Error: bedShear output was not as expected.")
		} else {
			err = os.Remove(v.OutFile)
			exception.PanicOnErr(err)
		}
	}
}
