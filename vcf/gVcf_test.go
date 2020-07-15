package vcf

import (
	"github.com/vertgenlab/gonomics/fileio"
	"testing"
)

func TestASFilter(t *testing.T) {
	file := fileio.EasyOpen("testdata/multiSampleTest.vcf")
	defer file.Close()
	header := ReadHeader(file)
	sampleHash := HeaderToMaps(header)

	reader := make(chan *Vcf)
	go ReadToChan(file, reader)

	var passFilter []*Vcf
	var parentalOne, parentalTwo, fOne int16 = sampleHash.GIndex["LITC"], sampleHash.GIndex["MATA"], sampleHash.GIndex["CL12_wgs_merged"]
	for each := range reader {
		if ASFilter(each, parentalOne, parentalTwo, fOne) {
			//PrintReOrder(each, []int16{parentalOne, parentalTwo, fOne})
			passFilter = append(passFilter, each)
		}
	}
	var currGt []GenomeSample
	for i := 0; i < len(passFilter); i++ {
		currGt = GetAlleleGenotype(passFilter[i])
		if IsHeterozygous(currGt[parentalOne]) || IsHeterozygous(currGt[parentalTwo]) || IsHomozygous(currGt[fOne]) {
			t.Errorf("Error: Parental genomes should be Homozygous and F1 should be Heterozygous...\n")
		}
	}
}

//TODO: Finish Sam split testing code
/*
func TestSamToAlleles(t *testing.T) {
	SnpSearch("testdata/CL12_wgsTest.sam", "testdata/multiSampleTest.vcf", "CL12_wgs_merged", "LITC", "MATA", "testdata/CL12_FINAL_test")
}*/
