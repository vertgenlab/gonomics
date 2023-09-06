package main

import (
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/motif"
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
	ResidualFilter     float64
}{
	{InFile: "testdata/STR012.fa",
		MatrixFile:         "testdata/jaspar.vertebrate.txt",
		ChromName:          "chr9",
		OutFile:            "testdata/tmp.tfMatchComp.bed",
		PropMatch:          0.8,
		MatrixFileType:     "Frequency",
		Pseudocounts:       0.1,
		RefStart:           113944,
		ExpectedFile:       "testdata/expected.tfMatchComp.bed",
		OutputAsProportion: true,
		ResidualFilter:     0.1,
	},
}

func TestTfMatchComp(t *testing.T) {
	var err error
	var s motif.MatchCompSettings
	for _, v := range TfMatchCompTests {
		s = motif.MatchCompSettings{
			MotifFile:          v.MatrixFile,
			ChromName:          v.ChromName,
			OutFile:            v.OutFile,
			PropMatch:          v.PropMatch,
			MotifType:          v.MatrixFileType,
			Pseudocounts:       v.Pseudocounts,
			RefStart:           v.RefStart,
			OutputAsProportion: v.OutputAsProportion,
			ResidualFilter:     v.ResidualFilter,
		}
		tfMatchComp(s, v.InFile)
		if !fileio.AreEqual(v.OutFile, v.ExpectedFile) {
			t.Errorf("Error in tfMatchComp. Output was not as expected.")
		} else {
			err = os.Remove(v.OutFile)
			exception.PanicOnErr(err)
		}
	}
}
