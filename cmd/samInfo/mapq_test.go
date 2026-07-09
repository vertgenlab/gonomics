package main

import (
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"os"
	"testing"
)

var MaqqTests = []struct {
	infile       string
	outhistogram string
	bedregions   string
	expected     string
}{
	{infile: "testdata/mapq/test1.sam", outhistogram: "testdata/mapq/out.hist.txt", bedregions: "", expected: "testdata/mapq/exp.hist.txt"},
	{infile: "testdata/mapq/test1.bam", outhistogram: "testdata/mapq/out.hist.txt", bedregions: "testdata/mapq/bed1.bed", expected: "testdata/mapq/exp.hist.txt"},
	{infile: "testdata/mapq/test1.bam", outhistogram: "testdata/mapq/out.hist.bed2.txt", bedregions: "testdata/mapq/bed2.bed", expected: "testdata/mapq/exp.hist.bed2.txt"},
}

func TestMapQ(t *testing.T) {
	var s mapqSettings

	for _, v := range MaqqTests {
		s = mapqSettings{
			InFile:     v.infile,
			OutFile:    v.outhistogram,
			BedRegions: v.bedregions,
		}
		mapq(s)
		if !fileio.AreEqual(v.expected, v.outhistogram) {
			t.Errorf("Error mapq: The exptected and observed files are not equal\n")
		} else {
			exception.PanicOnErr(os.Remove(v.outhistogram))
		}
	}
}
