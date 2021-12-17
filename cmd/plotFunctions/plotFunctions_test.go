package main

import (
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"os"
	"testing"
)

//note that the discrete functions are plotted in main, so these are tested in
//the popgen package.

var PlotFunctionsTests = []struct {
	Function     string
	FunctionArgs string
	Left         float64
	Right        float64
	Bins         int
	OutFile      string
	ExpectedFile string
}{
	{"AfsStationarity", "0.001", 0.001, 0.999, 100, "testdata/tmp.afsStationarity.txt", "testdata/expected.afsStationarity.txt"},
	{"Beta", "0.5,0.5", 0.001, 0.999, 100, "testdata/tmp.beta.txt", "testdata/expected.beta.txt"},
	{"Gamma", "0.5,0.5", 0.001, 0.999, 100, "testdata/tmp.gamma.txt", "testdata/expected.gamma.txt"},
	{"Normal", "0,0.5", -4, 4, 100, "testdata/tmp.normal.txt", "testdata/expected.normal.txt"},
}

func TestPlotFunctions(t *testing.T) {
	var err error
	for _, v := range PlotFunctionsTests {
		plotContinuousFunctions(v.Function, v.FunctionArgs, v.Left, v.Right, v.Bins, v.OutFile)
		if !fileio.AreEqual(v.OutFile, v.ExpectedFile) {
			t.Errorf("Error in plotFunctions, continuous functions. Output did not match expected.")
		} else {
			err = os.Remove(v.OutFile)
			exception.PanicOnErr(err)
		}
	}
}
