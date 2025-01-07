package main

import (
	"os"
	"testing"

	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
)

var GcContentTests = []struct {
	BedFile      string
	FaFile       string
	OutFile      string
	ExpectedFile string
	MultiFaMode  bool
	Species      string
}{
	{BedFile: "testdata/test.bed",
		FaFile:       "testdata/test.fa",
		OutFile:      "testdata/tmp.out.bed",
		ExpectedFile: "testdata/expected.bed",
		MultiFaMode:  false,
		Species:      ""},
	{BedFile: "testdata/multiFa.bed",
		FaFile:       "testdata/multiFa.fa",
		OutFile:      "testdata/tmp.multiFa.out.bed",
		ExpectedFile: "testdata/expected.multiFa.bed",
		MultiFaMode:  true,
		Species:      "Human_Chimp_Ancestor"},
}

func TestGcContent(t *testing.T) {
	var err error
	for _, v := range GcContentTests {
		gcContent(v.BedFile, v.FaFile, v.OutFile, v.MultiFaMode, v.Species)
		if !fileio.AreEqual(v.OutFile, v.ExpectedFile) {
			t.Errorf("Error in gcContent. Output did not match expected.")
		} else {
			err = os.Remove(v.OutFile)
			exception.PanicOnErr(err)
		}
	}
}
