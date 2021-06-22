package main

import (
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"os"
	"testing"
)

var VcfAfsTests = []struct {
	inputFile               string
	expectedFile            string
	unPolarized             bool
	plotSelectionLikelihood string
	leftBound               float64
	rightBound              float64
	numberOfPoints          int
	integralError           float64
}{
	{"testdata/simulate.N100.S100.Seed19.Alpha0.01.vcf", "testdata/test.afs.txt", true, "", -10, 10, 100, 1e-5},
}

func TestVcfAfs(t *testing.T) {
	var err error
	for _, v := range VcfAfsTests {
		s := Settings{
			UnPolarized:             v.unPolarized,
			PlotSelectionLikelihood: v.plotSelectionLikelihood,
			LeftBound:               v.leftBound,
			RightBound:              v.rightBound,
			NumberOfPoints:          v.numberOfPoints,
			IntegralError:           v.integralError,
		}
		vcfAfs(v.inputFile, "tmp.txt", s)
		if !fileio.AreEqual("tmp.txt", v.expectedFile) {
			t.Errorf("Error in VcfAfs.")
		} else {
			err = os.Remove("tmp.txt")
			exception.PanicOnErr(err)
		}
	}
}
