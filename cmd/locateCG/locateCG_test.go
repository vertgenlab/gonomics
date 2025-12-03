package main

import (
	"os"
	"testing"

	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
)

var LocateCGTest = []struct {
	inFa        string
	chromName   string
	outfile     string
	compare     bool
	cgtype      string
	expectedTxt string
}{
	{"testdata/single_genome_test.fa", "chr8", "chr8CG.bed", false, "", "testdata/single_genome_expected.bed"},
	{"testdata/twogenome_compare_test.fa", "chr8", "chr8gain.txt", true, "gain", "testdata/twogenome_gain_expected.txt"},
	{"testdata/twogenome_compare_test.fa", "chr8", "chr8loss.txt", true, "loss", "testdata/twogenome_loss_expected.txt"},
	{"testdata/twogenome_compare_test.fa", "chr8", "chr8cons.txt", true, "cons", "testdata/twogenome_cons_expected.txt"},
}

func TestLocateCG(t *testing.T) {
	var err error
	var s Settings
	for _, v := range LocateCGTest {
		s = Settings{
			InFa:      v.inFa,
			ChromName: v.chromName,
			Outfile:   v.outfile,
			Compare:   v.compare,
			CGtype:    v.cgtype,
		}
		if s.Compare {
			compareCG(s)
		} else {
			locateCG(s)
		}
		if !fileio.AreEqual(v.expectedTxt, v.outfile) {
			t.Errorf("Error in locateCG.go, expected output file != actual output file")
		} else {
			err = os.Remove(v.outfile)
			exception.PanicOnErr(err)
		}
	}
}
