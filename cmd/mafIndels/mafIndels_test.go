package main

import (
	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/exception"
	"os"
	"testing"
)

// use struct to specify testdata.
var mafIndelsTests = []struct {
	inMaf      string
	speciesIns string
	speciesDel string
	threshold  float64
	outInsBed  string
	outDelBed  string
}{
	{"testdata/in_hg38_vs_rheMac10_1.maf", "hg38", "rheMac10", 0.1, "testdata/outIns_hg38_1.bed", "testdata/outDel_rheMac10_1.bed"},
	{"testdata/in_hg38_vs_rheMac10_2.maf", "hg38", "rheMac10", 0.1, "testdata/outIns_hg38_2.bed", "testdata/outDel_rheMac10_2.bed"},
}

func TestMafIndels(t *testing.T) {
	for _, v := range mafIndelsTests {
		mafIndels(v.inMaf, v.speciesIns, v.speciesDel, v.threshold, "outIns_tmp.bed", "outDel_tmp.bed")

		records_ins := bed.Read("outIns_tmp.bed")
		expected_ins := bed.Read(v.outInsBed)
		if !bed.AllAreEqual(records_ins, expected_ins) {
			t.Errorf("Error: processing mafIndels for outIns.")
		}
		err := os.Remove("outIns_tmp.bed")
		exception.PanicOnErr(err)

		// because outDel might have illegal bed entries (start==end) that will cause bed.Read to log fatal, use fileio to compare test and expected
		if !fileio.AreEqual("outDel_tmp.bed", v.outDelBed) {
			t.Errorf("Error in mafIndels for outDel.")
		} else {
			err = os.Remove("outDel_tmp.bed")
			exception.PanicOnErr(err)
		}
	}
}
