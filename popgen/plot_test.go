package popgen

import (
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"os"
	"testing"
)

func BenchmarkAfsF1e7(b *testing.B) {
	for i := 0; i < b.N; i++ {
		PlotAfsF(0.03, 1322, "/dev/null", 1e-7)
	}
}

func BenchmarkAfsF1e6(b *testing.B) {
	for i := 0; i < b.N; i++ {
		PlotAfsF(0.03, 1322, "/dev/null", 1e-6)
	}
}

func BenchmarkAfsF1e5(b *testing.B) {
	for i := 0; i < b.N; i++ {
		PlotAfsF(0.03, 1322, "/dev/null", 1e-5)
	}
}

func BenchmarkAfsF1e4(b *testing.B) {
	for i := 0; i < b.N; i++ {
		PlotAfsF(0.03, 1322, "/dev/null", 1e-4)
	}
}

var PlotAfsFTests = []struct {
	alpha         float64
	n             int
	outFile       string
	integralError float64
	expectedFile  string
}{
	{0.01, 10, "testdata/tmp.AfsF.txt", 1e-5, "testdata/expected.AfsF.txt"},
}

func TestPlotAfsF(t *testing.T) {
	var err error
	for _, v := range PlotAfsFTests {
		PlotAfsF(v.alpha, v.n, v.outFile, v.integralError)
		if !fileio.AreEqual(v.outFile, v.expectedFile) {
			t.Errorf("Error in PlotAfsF. Output did not match expected.")
		} else {
			err = os.Remove(v.outFile)
			exception.PanicOnErr(err)
		}
	}
}

var PlotAfsPmfTests = []struct {
	alpha         float64
	n             int
	outFile       string
	integralError float64
	derived       bool
	ancestral     bool
	expectedFile  string
}{
	{0.01, 10, "testdata/tmp.AfsPmf.txt", 1e-5, false, false, "testdata/expected.AfsPmf.txt"},
	{0.01, 10, "testdata/tmp.AfsPmfDerived.txt", 1e-5, true, false, "testdata/expected.AfsPmfDerived.txt"},
	{0.01, 10, "testdata/tmp.AfsPmfAncestral.txt", 1e-5, false, true, "testdata/expected.AfsPmfAncestral.txt"},
}

func TestPlotAfsPmf(t *testing.T) {
	var err error
	for _, v := range PlotAfsPmfTests {
		PlotAfsPmf(v.alpha, v.n, v.outFile, v.integralError, v.derived, v.ancestral)
		if !fileio.AreEqual(v.outFile, v.expectedFile) {
			t.Errorf("Error in PlotAfsPmf. Output was not as expected.")
		} else {
			err = os.Remove(v.outFile)
			exception.PanicOnErr(err)
		}
	}
}

var PlotAfsLikelihoodTests = []struct {
	vcfFile                 string
	outFile                 string
	leftBound               float64
	rightBound              float64
	numPoints               int
	integralError           float64
	divergenceAscertainment bool
	D                       int
	expectedFile            string
}{
	{"testdata/simulated.alpha4.N100.S100.seed19.vcf",
		"testdata/tmp.likelihood.txt",
		-9,
		9,
		21,
		1e-5,
		false,
		1,
		"testdata/expected.likelihoodPlot.txt",
	},
}

func TestPlotAfsLikelihood(t *testing.T) {
	var err error
	var currAfs *Afs
	for _, v := range PlotAfsLikelihoodTests {
		currAfs, err = VcfToAfs(v.vcfFile, false, false, false)
		exception.PanicOnErr(err)
		PlotAfsLikelihood(*currAfs, v.outFile, v.leftBound, v.rightBound, v.numPoints, v.integralError, v.divergenceAscertainment, v.D)
		if !fileio.AreEqual(v.outFile, v.expectedFile) {
			t.Errorf("Error in PlotAfsLikelihood. Output not as expected.")
		} else {
			err = os.Remove(v.outFile)
			exception.PanicOnErr(err)
		}
	}
}

var PlotAscertainmentProbabilityTests = []struct {
	outFile               string
	n                     int
	d                     int
	expectedDerivedFile   string
	expectedAncestralFile string
}{
	{"testdata/tmp.AscertainmentProbability.txt",
		50,
		1,
		"testdata/expected.DerivedAscertainmentProbability.txt",
		"testdata/expected.AncestralAscertainmentProbability.txt",
	},
}

func TestPlotAscertainmentProbabilityNumerators(t *testing.T) {
	var err error
	for _, v := range PlotAscertainmentProbabilityTests {
		PlotDerivedAscertainmentProbability(v.outFile, v.n, v.d)
		if !fileio.AreEqual(v.outFile, v.expectedDerivedFile) {
			t.Errorf("Error in PlotDerivedAscertainmentProbability. Output not as expected.")
		} else {
			err = os.Remove(v.outFile)
			exception.PanicOnErr(err)
		}

		PlotAncestralAscertainmentProbability(v.outFile, v.n, v.d)
		if !fileio.AreEqual(v.outFile, v.expectedAncestralFile) {
			t.Errorf("Error in PlotAncestralAscertainmentProbability. Output not as expected.")
		} else {
			err = os.Remove(v.outFile)
			exception.PanicOnErr(err)
		}
	}
}

var PlotAscertainmentDenominatorProbabilitiesTests = []struct {
	outFile           string
	n                 int
	d                 int
	alpha             float64
	integralError     float64
	expectedDerived   string
	expectedAncestral string
}{
	{"testdata/tmp.Denominator.txt",
		10,
		1,
		0.01,
		1e-5,
		"testdata/expected.DerivedDenominator.txt",
		"testdata/expected.AncestralDenominator.txt",
	},
}

func TestPlotAscertainmentDenominator(t *testing.T) {
	var err error
	for _, v := range PlotAscertainmentDenominatorProbabilitiesTests {
		PlotDerivedAscertainmentDenominator(v.outFile, v.n, v.d, v.alpha, v.integralError)
		if !fileio.AreEqual(v.outFile, v.expectedDerived) {
			t.Errorf("Error in PlotDerivedAscertainmentDenominator.")
		} else {
			err = os.Remove(v.outFile)
			exception.PanicOnErr(err)
		}

		PlotAncestralAscertainmentDenominator(v.outFile, v.n, v.d, v.alpha, v.integralError)
		if !fileio.AreEqual(v.outFile, v.expectedAncestral) {
			t.Errorf("Error in PlotAncestralAscertainmentDenominator.")
		} else {
			err = os.Remove(v.outFile)
			exception.PanicOnErr(err)
		}
	}
}
