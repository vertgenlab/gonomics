package main

import (
	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/common"
	"os"
	"testing"
)

//use struct to specify testdata
var mafIndelsTests = []struct {
	in_maf      string
	species_ins string
	species_del string
	threshold   float64
	outIns_bed  string
	outDel_bed  string
}{
	{"testdata/in_hg38_vs_rheMac10_1.maf", "hg38", "rheMac10", 0.1, "testdata/outIns_hg38_1.bed", "testdata/outDel_rheMac10_1.bed"},
	{"testdata/in_hg38_vs_rheMac10_2.maf", "hg38", "rheMac10", 0.1, "testdata/outIns_hg38_2.bed", "testdata/outDel_rheMac10_2.bed"},
}

func TestMafIndels(t *testing.T) {
	for _, v := range mafIndelsTests {
		mafIndels(v.in_maf, v.species_ins, v.species_del, v.threshold, "outIns_tmp.bed", "outDel_tmp.bed")

		records_ins := bed.Read("outIns_tmp.bed")
		expected_ins := bed.Read(v.outIns_bed)
		if !bed.AllAreEqual(records_ins, expected_ins) {
			t.Errorf("Error in mafIndels for outIns.")
		}
		err := os.Remove("outIns_tmp.bed")
		if err != nil {
			common.ExitIfError(err)
		}

		records_del := bed.Read("outDel_tmp.bed")
		expected_del := bed.Read(v.outDel_bed)
		if !bed.AllAreEqual(records_del, expected_del) {
			t.Errorf("Error in mafIndels for outDel.")
		}
		err = os.Remove("outDel_tmp.bed")
		if err != nil {
			common.ExitIfError(err)
		}

	}
}
