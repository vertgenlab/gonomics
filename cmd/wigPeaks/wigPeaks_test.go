package main

import (
	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/common"
	"os"
	"testing"
)

//use struct to specify testdata
var wigPeaksTests = []struct {
	in_wig    string
	out_bed   string
	threshold float64
}{
	{"testdata/in_wig_1.wig", "testdata/out_bed_1.bed", 20},
	{"testdata/in_wig_2.wig", "testdata/out_bed_2.bed", 50},
}

func TestWigPeaks(t *testing.T) {
	for _, v := range wigPeaksTests {
		wigPeaks(v.in_wig, "out_bed_tmp.bed", v.threshold)

		records := bed.Read("out_bed_tmp.bed")
		expected := bed.Read(v.out_bed)
		if !bed.AllAreEqual(records, expected) {
			t.Errorf("Error in wigPeaks.")
		}
		//DEBUG: don't delete out_bed_tmp.bed, unimport common, os
		err := os.Remove("out_bed_tmp.bed")
		if err != nil {
			common.ExitIfError(err)
		}
	}
}
