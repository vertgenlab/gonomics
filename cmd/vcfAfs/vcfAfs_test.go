package main

import (
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"os"
	"testing"
)

var VcfAfsTests = []struct {
	inputFile               string
	outputFile string
	expectedFile            string
	unPolarized             bool
	plotSelectionLikelihood string
	leftBound               float64
	rightBound              float64
	numberOfPoints          int
	integralError           float64
	includeRef	bool
}{
	{"testdata/simulate.N100.S100.Seed19.Alpha0.01.vcf",
		"testdata/tmp.txt",
		"testdata/expected.afs.txt",
		false,
		"",
		-10,
		10,
		100,
		1e-5,
		false},
	{"testdata/simulate.N100.S100.Seed19.Alpha0.01.vcf",
		"testdata/tmp.IncludeRef.txt",
		"testdata/expected.IncludeRef.afs.txt",
		false,
		"",
		-10,
		10,
		100,
		1e-5,
		true},
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
			IncludeRef: v.includeRef,
		}
		vcfAfs(v.inputFile, v.outputFile, s)
		if !fileio.AreEqual(v.outputFile, v.expectedFile) {
			t.Errorf("Error in VcfAfs.")
		} else {
			err = os.Remove(v.outputFile)
			exception.PanicOnErr(err)
		}
	}
}
