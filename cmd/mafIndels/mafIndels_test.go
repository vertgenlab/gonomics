package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/maf"
	"log"
	"testing"
)

//use struct to specify testdata
var mafIndelsTests = []struct {
	in_maf          string
	species_ins 		string
	species_del     string
	outIns_bed      string
	outDel_bed			string
}{
	{"testdata/in_hg38_vs_rheMac10.maf", "hg38", "rheMac10", "testdata/outIns_hg38.bed", "testdata/outDel_rheMac10.bed"},
	{"testdata/in_hg38_vs_rheMac10_2.maf", "hg38", "rheMac10", "testdata/outIns_hg38_2.bed", "testdata/outDel_rheMac10_2.bed"},
}

func TestMafIndels(t *testing.T) {
	for _, v := range mafIndelsTests {

		//TODO: edit code from here
		if vcf.IsVcfFile(v.inputFile) {
			lift(v.chainFile, v.inputFile, "tmp.vcf", v.faFile, "tmp.unmapped", 0.95)
			records := vcf.Read("tmp.vcf")
			expected := vcf.Read(v.expectedOutputFile)
			if !vcf.AllEqual(records, expected) {
				t.Errorf("Error in Lift for vcf.")
			}
			err := os.Remove("tmp.vcf")
			if err != nil {
				common.ExitIfError(err)
			}
			err = os.Remove("tmp.unmapped")
			if err != nil {
				common.ExitIfError(err)
			}
		} else {
			lift(v.chainFile, v.inputFile, "tmp.bed", v.faFile, "tmp.unmapped", 0.95)
			records := bed.Read("tmp.bed")
			expected := bed.Read(v.expectedOutputFile)
			if !bed.AllAreEqual(records, expected) {
				t.Errorf("Error in Lift for bed.")
			}
			err := os.Remove("tmp.bed")
			if err != nil {
				common.ExitIfError(err)
			}
			err = os.Remove("tmp.unmapped")
			if err != nil {
				common.ExitIfError(err)
			}
		}
	}
}
