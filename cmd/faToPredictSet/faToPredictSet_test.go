package main

import (
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"os"
	"testing"
)

var FaToPredictSetTests = []struct {
	InFile       string
	OutFile      string
	ExpectedFile string
	WindowSize   int
	Stride       int
	WithRevComp  bool
}{
	{InFile: "testdata/test.fa",
		OutFile:      "testdata/tmp.txt",
		ExpectedFile: "testdata/expected.txt",
		WindowSize:   10,
		Stride:       1,
		WithRevComp:  false,
	},
	{InFile: "testdata/test.fa",
		OutFile:      "testdata/tmp.withRevComp.txt",
		ExpectedFile: "testdata/expected.withRevComp.txt",
		WindowSize:   10,
		Stride:       1,
		WithRevComp:  true,
	},
}

func TestFaToPredictSet(t *testing.T) {
	var err error
	for _, v := range FaToPredictSetTests {
		s := Settings{
			InFile:      v.InFile,
			OutFile:     v.OutFile,
			WindowSize:  v.WindowSize,
			Stride:      v.Stride,
			WithRevComp: v.WithRevComp,
		}
		faToPredictSet(s)
		if !fileio.AreEqual(v.ExpectedFile, v.OutFile) {
			t.Errorf("Error in faToPredictSet. Output was not as expected.")
		} else {
			err = os.Remove(v.OutFile)
			exception.PanicOnErr(err)
		}
	}
}
