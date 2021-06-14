package main

import (
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"os"
	"testing"
)

var McmcTraceStatsTests = []struct {
	InFile        string
	OutFile       string
	ExpectedFile  string
	HdiProportion float64
	BurnIn        int
	ParameterName string
}{
	{"testdata/Rand.trace.txt", "tmp.txt", "testdata/Rand.trace.stats.txt", 0.95, 5000, "Mu"},
	{"testdata/Rand.trace.txt", "tmp.txt", "testdata/Rand.trace.stats.sigma.txt", 0.95, 5000, "Sigma"},
}

func TestMcmcTraceStats(t *testing.T) {
	var err error
	for _, v := range McmcTraceStatsTests {
		mcmcTraceStats(v.InFile, v.OutFile, v.HdiProportion, v.BurnIn, v.ParameterName)
		if !fileio.AreEqual(v.OutFile, v.ExpectedFile) {
			t.Errorf("Error in mcmcTraceStats.")
		}
		err = os.Remove(v.OutFile)
		exception.PanicOnErr(err)
	}
}
