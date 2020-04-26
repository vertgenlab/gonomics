package vcf

import (
	"fmt"
	"testing"
)

func TestVcfToGvcf(t *testing.T) {
	dict := HeaderToMaps("testdata/small.LITCxMATA.GENOTYPED.vcf")
	vcfPipe := make(chan *Vcf)
	go ReadToChan("testdata/small.LITCxMATA.GENOTYPED.vcf", vcfPipe)
	for record := range vcfPipe {
		if SimpleASFilter(record, dict.HapIdx["MATAxLITC_cl12w16-7_L6"], dict.HapIdx["LITC"], dict.HapIdx["MATA"]) {
			ViewGenotypeVcf(record)
		}
	}
}

func TestSampleSheetFilter(t *testing.T) {
	dict := HeaderToMaps("testdata/small.LITCxMATA.GENOTYPED.vcf")
	vcfPipe := make(chan *Vcf)
	Aa, aa := ReadFilterList("testdata/sampleSheet.csv", dict.HapIdx)
	index1, index2 := MapNameToIndex(dict.HapIdx, Aa), MapNameToIndex(dict.HapIdx, aa)
	go ReadToChan("testdata/small.LITCxMATA.GENOTYPED.vcf", vcfPipe)
	fmt.Printf("#CHROM\tPOS\tREF\tALT\tCL12merged.bam\tMATAxLITC_cl12w16-3_L6\tMATAxLITC_cl12w16-7_L6\tLITC\tMATA\n")
	for record := range vcfPipe {
		if ASFilter(record, index1, index2) {
			prettyShuffle(record, index1, index2)
		}
	}
}

func TestSamToAlleles(t *testing.T) {
	SnpSearch("testdata/CL12.wgs.merged.scaffold_1_30030-88166.sam", "testdata/LITCxMATA.trio.test.vcf", "testdata/sampleSheet.csv")
}
