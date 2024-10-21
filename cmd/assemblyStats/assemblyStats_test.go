package main

import (
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"os"
	"testing"
)

var assemblyStatsTests = []struct {
	inFa            string
	outFile         string
	expFile         string
	lowercaseAsGaps bool
}{
	{"testdata/test.fa", "testdata/outFalse.txt", "testdata/expFalse.txt", false},
	{"testdata/test.fa", "testdata/outTrue.txt", "testdata/expTrue.txt", true},
}

func TestAssemblyStats(t *testing.T) {
	var err error
	for _, v := range assemblyStatsTests {
		assemblyStats(v.inFa, v.outFile, v.lowercaseAsGaps)
		if !fileio.AreEqual(v.outFile, v.expFile) {
			t.Errorf("Error in assemblyStats")
		} else {
			err = os.Remove(v.outFile)
			exception.PanicOnErr(err)
		}
	}
}
