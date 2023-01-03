package main

import (
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"os"
	"testing"
)

var TfMatchCompTests = []struct {
	InFile             string
	MatrixFile         string
	ChromName          string
	OutFile            string
	PropMatch          float64
	MatrixFileType     string
	Pseudocounts       float64
	RefStart           int
	ExpectedFile       string
	OutputAsProportion bool
}{
	{InFile: "testdata/STR012.fa",
		MatrixFile:         "testdata/jaspar.vertebrate.pfm",
		ChromName:          "chr9",
		OutFile:            "testdata/tmp.tfMatchComp.bed",
		PropMatch:          0.8,
		MatrixFileType:     "Frequency",
		Pseudocounts:       0.1,
		RefStart:           113944,
		ExpectedFile:       "testdata/expected.tfMatchComp.bed",
		OutputAsProportion: false},
	{InFile: "testdata/STR012.fa",
		MatrixFile:         "testdata/jaspar.vertebrate.pfm",
		ChromName:          "chr9",
		OutFile:            "testdata/tmp.outputProp.tfMatchComp.bed",
		PropMatch:          0.8,
		MatrixFileType:     "Frequency",
		Pseudocounts:       0.1,
		RefStart:           113944,
		ExpectedFile:       "testdata/expected.outputProp.tfMatchComp.bed",
		OutputAsProportion: true},
}

func TestTfMatchComp(t *testing.T) {
	var err error
	var s Settings
	for _, v := range TfMatchCompTests {
		s = Settings{
			InFile:         v.InFile,
			MatrixFile:     v.MatrixFile,
			ChromName:      v.ChromName,
			OutFile:        v.OutFile,
			PropMatch:      v.PropMatch,
			MatrixFileType: v.MatrixFileType,
			Pseudocounts:   v.Pseudocounts,
			RefStart:       v.RefStart,
			OutputAsProportion: v.OutputAsProportion,
		}
		tfMatchComp(s)
		if !fileio.AreEqual(v.OutFile, v.ExpectedFile) {
			t.Errorf("Error in tfMatchComp. Output was not as expected.")
		} else {
			err = os.Remove(v.OutFile)
			exception.PanicOnErr(err)
		}
	}
}
