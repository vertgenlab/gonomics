package main

import (
	"os"
	"testing"

	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
)

var CountPairOfBasesTest = []struct {
	faDir       string
	bed         string
	baseone     string
	basetwo     string
	outfile     string
	compare     bool
	expectedTxt string
}{
	{"./testdata/singlegenome", "./testdata/single_genome_test.bed", "C", "G", "single_genome_out.bed", false, "./testdata/single_genome_expected.bed"},
	{"./testdata/twogenome", "./testdata/two_genome_test.bed", "C", "G", "two_genome_out.bed", true, "./testdata/two_genome_expected.txt"},
}

func TestCountPairOfBases(t *testing.T) {
	var err error
	var s Settings
	for _, v := range CountPairOfBasesTest {
		s = Settings{
			FaDir:   v.faDir,
			Bed:     v.bed,
			BaseOne: v.baseone,
			BaseTwo: v.basetwo,
			Outfile: v.outfile,
			Compare: v.compare,
		}
		countByChrom(s)

		if !fileio.AreEqual(v.expectedTxt, v.outfile) {
			t.Errorf("Error in countPairOfBases.go, expected output file != actual output file")
		} else {
			err = os.Remove(v.outfile)
			exception.PanicOnErr(err)
		}
	}
}
