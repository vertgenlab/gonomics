package vcf

import (
	"testing"
)

func TestASFilter(t *testing.T) {
	reader, header := GoReadToChan("testdata/multiSampleTest.vcf")
	sampleHash := HeaderToMaps(header)

	var passFilter []Vcf
	var parentalOne, parentalTwo, fOne int16 = sampleHash.GIndex["LITC"], sampleHash.GIndex["MATA"], sampleHash.GIndex["CL12_wgs_merged"]
	for each := range reader {
		if ASFilter(each, parentalOne, parentalTwo, fOne) {
			//PrintReOrder(each, []int16{parentalOne, parentalTwo, fOne})
			passFilter = append(passFilter, each)
		}
	}
	var currGt []GenomeSample
	for i := 0; i < len(passFilter); i++ {
		currGt = passFilter[i].Samples
		if IsHeterozygous(currGt[parentalOne]) || IsHeterozygous(currGt[parentalTwo]) || IsHomozygous(currGt[fOne]) {
			t.Errorf("Error: Parental genomes should be Homozygous and F1 should be Heterozygous...\n")
		}
	}
}
