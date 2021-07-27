package main

import (
	"testing"
	"os"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/exception"
)

var bedToWigTests = []struct {
	inFile string
	refFile string
	outFile string
	expectedFile string
	method string
	missing float64
}{
	{"testdata/test.bed", "testdata/ref.chrom.sizes", "testdata/test.Score.wig", "testdata/score.Expected.wig", "Score", 0},
	{"testdata/test.bed", "testdata/ref.chrom.sizes", "testdata/test.Reads.wig", "testdata/reads.Expected.wig", "Reads", 0},
	{"testdata/test.bed", "testdata/ref.chrom.sizes", "testdata/test.Name.wig", "testdata/name.Expected.wig", "Name", 0},
	{"testdata/test.bed", "testdata/ref.chrom.sizes", "testdata/test.missing.Name.wig", "testdata/name.missing.Expected.wig", "Name", -1.0},
}

func TestBedToWig(t *testing.T) {
	var err error
	for _, v := range bedToWigTests {
		bedToWig(v.method, v.inFile, v.refFile, v.outFile, v.missing)
		if !fileio.AreEqual(v.outFile, v.expectedFile) {
			t.Errorf("Error in BedToWig.")
		} else {
			err = os.Remove(v.outFile)
			exception.PanicOnErr(err)
		}
	}
}
