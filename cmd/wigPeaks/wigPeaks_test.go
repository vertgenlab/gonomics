package main

import (
	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/exception"
	"os"
	"testing"
)

//use struct to specify testdata
var wigPeaksTests = []struct {
	inWig       string
	outBed      string
	expectedBed string
	threshold   float64
	findMinima  bool
}{
	//{"testdata/in_wig_1.wig", "testdata/out_bed_tmp.bed", "testdata/out_bed_1.bed", 20, false},
	//{"testdata/in_wig_2.wig", "testdata/out_bed_tmp.bed", "testdata/out_bed_2.bed", 50, false},
	{"testdata/in_wig_1.wig", "testdata/tmp.Minima.bed", "testdata/expected.minima.bed", 50, true},
}

func TestWigPeaks(t *testing.T) {
	var err error
	for _, v := range wigPeaksTests {
		wigPeaks(v.inWig, v.outBed, v.threshold, v.findMinima)
		records := bed.Read(v.outBed)
		expected := bed.Read(v.expectedBed)
		if !bed.AllAreEqual(records, expected) {
			t.Errorf("Error in wigPeaks, created bed and test bed are not equal.")
		} else {
			err = os.Remove(v.outBed)
			exception.PanicOnErr(err)
		}
	}
}
