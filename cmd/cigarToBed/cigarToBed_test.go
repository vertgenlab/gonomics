package main

import (
	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/common"
	"github.com/vertgenlab/gonomics/fileio"
	"os"
	"testing"
)

var CigarToBedTests = []struct {
	inputFileOne_Name   string
	inputFileTwo_Name   string
	outFa               string
	FirstPos_InsBed     int
	FirstPos_DelBed     int
	Chrom               string
	outIns_bed_expected string
	outDel_bed_expected string
}{
	{"testdata/sethvsraven/seth.fa", "testdata/sethvsraven/raven.fa", "", 1, 1, "chr1", "testdata/sethvsraven/affineGap_sethvsraven_ins.bed", "testdata/sethvsraven/affineGap_sethvsraven_del.bed"},
	{"testdata/firstTest/testRegion10kb_PanTro6.fa", "testdata/firstTest/testRegion10kb_hg38.fa", "", 119320000, 116703287, "chr1", "testdata/firstTest/affineGap_PanTro6vshg38_ins.bed", "testdata/firstTest/affineGap_PanTro6vshg38_del.bed"},
}

func TestCigarToBed(t *testing.T) {
	for _, v := range CigarToBedTests {
		inputFileOne := fileio.EasyOpen(v.inputFileOne_Name)
		inputFileTwo := fileio.EasyOpen(v.inputFileTwo_Name)
		GlobalAlignment_CigarToBed(inputFileOne, inputFileTwo, v.outFa, "ins_tmp.bed", "del_tmp.bed", v.FirstPos_InsBed, v.FirstPos_DelBed, v.Chrom)

		records_ins := bed.Read("ins_tmp.bed")
		expected_ins := bed.Read(v.outIns_bed_expected)
		if !bed.AllAreEqual(records_ins, expected_ins) {
			t.Errorf("Error in mafIndels for outIns.")
		}
		err := os.Remove("ins_tmp.bed")
		if err != nil {
			common.ExitIfError(err)
		}

		records_del := bed.Read("del_tmp.bed")
		expected_del := bed.Read(v.outDel_bed_expected)
		if !bed.AllAreEqual(records_del, expected_del) {
			t.Errorf("Error in mafIndels for outDel.")
		}
		err = os.Remove("del_tmp.bed")
		if err != nil {
			common.ExitIfError(err)
		}

	}
}
