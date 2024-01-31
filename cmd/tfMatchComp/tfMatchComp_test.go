package main

import (
	"github.com/vertgenlab/gonomics/exception"
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
	GcContent          float64
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
		GcContent:          0.5,
	},
}

func TestTfMatchComp(t *testing.T) {
	var err error
	var s motif.MatchCompSettings
	var tolerance float64 = 0.1
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
			GcContent:          v.GcContent,
		}
		tfMatchComp(s, v.InFile)
		if !motif.AlmostEqualTest(v.ExpectedFile, v.OutFile, tolerance) {
			t.Errorf("Error: Motif are not equal within a tolorance %v...", tolerance)
		} else {
			err = os.Remove(v.OutFile)
			exception.PanicOnErr(err)
		}
	}
}
