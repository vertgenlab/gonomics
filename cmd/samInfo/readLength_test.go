package main

import (
	"github.com/vertgenlab/gonomics/fileio"
	"testing"
)

var ReadLengthTests = []struct {
	InFile       string
	OutFile      string
	ExpectedFile string
}{
	{InFile: "testdata/readLength/small.sam",
		OutFile:      "testdata/readLength/test.readLength.txt",
		ExpectedFile: "testdata/readLength/expected.readLength.txt",
	},
}

func TestReadLength(t *testing.T) {
	var s ReadLengthSettings
	for _, v := range ReadLengthTests {
		s = ReadLengthSettings{
			InFile:  v.InFile,
			OutFile: v.OutFile,
		}
		readLength(s)
		if !fileio.AreEqual(v.OutFile, v.ExpectedFile) {
			t.Errorf("Error: samInfo readLength output is not as expected.\n")
		} else {
			fileio.MustRemove(v.OutFile)
		}
	}
}
