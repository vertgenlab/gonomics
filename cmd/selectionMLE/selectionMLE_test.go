package main

import (
	//"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/popgen"
	"testing"
)

var SelectionMleTests = []struct {
	inputFile string
	outFile string
	expectedOutputFile string
	left float64
	right float64
	error float64
	unPolarized bool
	divergenceAscertainment bool
	integralError float64
	verbose int
}{
	{"testdata/simulated.alpha4.N100.S100.seed19.vcf", "tmp.txt", "testdata/expected4.txt", -10, 10, 1e-5, true, false, 1e-5, 0},
}

func TestSelectionMle(t *testing.T) {
	//var err error
	for _, v := range SelectionMleTests {
		s := popgen.MleSettings{
			Left: v.left,
			Right: v.right,
			Error: v.error,
			UnPolarized: v.unPolarized,
			DivergenceAscertainment: v.divergenceAscertainment,
			D: 1,
			IntegralError: v.integralError,
			Verbose: v.verbose,
		}
		selectionMLE(v.inputFile, v.outFile, s)
		if !fileio.AreEqual(v.expectedOutputFile, v.outFile) {
			t.Errorf("Error in SelectionMLE.")
		}
		//err = os.Remove(v.outFile)
		//exception.PanicOnErr(err)
	}
}