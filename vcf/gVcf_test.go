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
	var parentalOne, parentalTwo, fOne int16 = sampleHash.IndexAllele["LITC"], sampleHash.IndexAllele["MATA"], sampleHash.IndexAllele["CL12_wgs_merged"]
	for each := range reader {
		if ASFilter(each, parentalOne, parentalTwo, fOne) {
			//PrintReOrder(each, []int16{parentalOne, parentalTwo, fOne})
			passFilter = append(passFilter, each)
		}
	}
	var currGt []Sample
	for i := 0; i < len(passFilter); i++ {
		currGt = GetAlleleGenotype(passFilter[i])
		if isHeterozygous(currGt[parentalOne]) || isHeterozygous(currGt[parentalTwo]) || isHomozygous(currGt[fOne]) {
			t.Errorf("Error: Parental genomes should be Homozygous and F1 should be Heterozygous...\n")
		}
	}
}