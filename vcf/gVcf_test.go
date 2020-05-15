package vcf

import (
	"github.com/vertgenlab/gonomics/fileio"
	"testing"
)

func TestVcfToGvcf(t *testing.T) {
	file := fileio.EasyOpen("testdata/testGenotype.vcf")
	defer file.Close()
	dict := HeaderToMaps(file)
	vcfPipe := make(chan *Vcf)
	go ReadToChan(file, vcfPipe)
	for record := range vcfPipe {
		if SimpleASFilter(record, dict.HapIdx["MATAxLITC_cl12w16-7_L6"], dict.HapIdx["LITC"], dict.HapIdx["MATA"]) {
			ViewGenotypeVcf(record)
		}
	}
}

func TestSamToAlleles(t *testing.T) {
	SnpSearch("testdata/CL12.wgs.merged.scaffold_1_30030-88166.sam", "testdata/LITCxMATA.trio.test.vcf", "CL12merged.bam", "LITC", "MATA", "testdata/CL12.wgs.merged")
}
