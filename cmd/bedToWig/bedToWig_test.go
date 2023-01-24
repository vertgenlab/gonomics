package main

import (
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"os"
	"testing"
)

var bedToWigTests = []struct {
	inFile       string
	refFile      string
	outFile      string
	expectedFile string
	method       string
	missing      float64
	useRange bool
}{
	{inFile: "testdata/test.bed",
		refFile: "testdata/ref.chrom.sizes",
		outFile: "testdata/test.Score.wig",
		expectedFile: "testdata/score.Expected.wig",
		method: "Score",
		missing: 0,
		useRange: false},
	{inFile: "testdata/test.bed",
		refFile: "testdata/ref.chrom.sizes",
		outFile: "testdata/test.Reads.wig",
		expectedFile: "testdata/reads.Expected.wig",
		method: "Reads",
		missing: 0,
		useRange: false},
	{inFile: "testdata/test.bed",
		refFile: "testdata/ref.chrom.sizes",
		outFile: "testdata/test.Name.wig",
		expectedFile: "testdata/name.Expected.wig",
		method: "Name",
		missing: 0,
		useRange: false},
	{inFile: "testdata/test.bed",
		refFile: "testdata/ref.chrom.sizes",
		outFile: "testdata/test.missing.Name.wig",
		expectedFile: "testdata/name.missing.Expected.wig",
		method: "Name",
		missing: -1.0,
		useRange: false},
	{inFile: "testdata/test.range.bed",
		refFile: "testdata/ref.chrom.sizes",
		outFile: "testdata/test.range.Name.wig",
		expectedFile: "testdata/name.range.Expected.wig",
		method: "Name",
		missing: -1.0,
		useRange: true},
	{inFile: "testdata/test.range.bed",
		refFile: "testdata/ref.chrom.sizes",
		outFile: "testdata/test.range.Score.wig",
		expectedFile: "testdata/score.range.Expected.wig",
		method: "Score",
		missing: -1.0,
		useRange: true},
}

func TestBedToWig(t *testing.T) {
	var err error
	for _, v := range bedToWigTests {
		bedToWig(v.method, v.inFile, v.refFile, v.outFile, v.missing, v.useRange)
		if !fileio.AreEqual(v.outFile, v.expectedFile) {
			t.Errorf("Error in BedToWig.")
		} else {
			err = os.Remove(v.outFile)
			exception.PanicOnErr(err)
		}
	}
}
