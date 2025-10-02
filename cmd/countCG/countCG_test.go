package main

import (
	"os"
	"testing"

	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
)

var CountCGTest = []struct {
	faDir       string
	bed         string
	outfile     string
	compare     bool
	expectedTxt string
}{
	{"./testdata/singlegenome", "./testdata/single_genome_test.bed", "single_genome_out.bed", false, "./testdata/single_genome_expected.bed"},
	{"./testdata/twogenome", "./testdata/two_genome_test.bed", "two_genome_out.bed", true, "./testdata/two_genome_expected.txt"},
}

func TestCountCG(t *testing.T) {
	var err error
	var s Settings
	for _, v := range CountCGTest {
		s = Settings{
			FaDir:   v.faDir,
			Bed:     v.bed,
			Outfile: v.outfile,
			Compare: v.compare,
		}
		countBychrom(s)

		if !fileio.AreEqual(v.expectedTxt, v.outfile) {
			t.Errorf("Error in countCG.go, expected output file != actual output file")
		} else {
			err = os.Remove(v.outfile)
			exception.PanicOnErr(err)
		}
	}
}
