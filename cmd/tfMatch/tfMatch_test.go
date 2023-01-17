package main

import (
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"os"
	"testing"
)

var tfMatchTests = []struct {
	InFile             string
	MatrixFile         string
	OutFile            string
	MatrixFileType     string
	PropMatch          float64
	Pseudocounts       float64
	OutputAsProportion bool
	ExpectedFile       string
}{
	{InFile: "testdata/chr1.upper.firstMb.fa",
		MatrixFile:         "testdata/jaspar.small.txt",
		OutFile:            "testdata/tmp.tfMatch.bed",
		MatrixFileType:     "Frequency",
		PropMatch:          0.8,
		Pseudocounts:       0.1,
		OutputAsProportion: false,
		ExpectedFile:       "testdata/expected.tfMatch.bed"},
}

func TestTfMatch(t *testing.T) {
	var err error
	var s Settings
	for _, v := range tfMatchTests {
		s = Settings{
			InFile:             v.InFile,
			MatrixFile:         v.MatrixFile,
			OutFile:            v.OutFile,
			MatrixFileType:     v.MatrixFileType,
			PropMatch:          v.PropMatch,
			Pseudocounts:       v.Pseudocounts,
			OutputasProportion: v.OutputAsProportion,
		}
		tfMatch(s)
		if !fileio.AreEqual(v.OutFile, v.ExpectedFile) {
			t.Errorf("Error in tfMatch. Output was not as expected.")
		} else {
			err = os.Remove(v.OutFile)
			exception.PanicOnErr(err)
		}
	}
}
