package alleles

import (
	"testing"
)

/*
func TestCountAlleles(t *testing.T) {
	samplemap := CountAlleles("testdata/human_chrM.fasta",
"testdata/human_chrM.sam", 0)
	if samplemap == nil {
		t.Errorf("Problem with CountAlleles")
	}
}

func TestGoCountAlleles(t *testing.T) {
	samplemap := GoCountAlleles("testdata/human_chrM.fasta",
"testdata/human_chrM.sam", 0, 10)
	if samplemap == nil {
		t.Errorf("Problem with GoCountAlleles")
	}
}

func TestAllelesToVcf(t *testing.T) {
	samplemap := CountAlleles("testdata/human_chrM.fasta",
"testdata/human_chrM.sam", 0)
	answer := AllelesToVcf(samplemap)
	if answer == nil {
		t.Errorf("Problem with AllelesToVcf")
	}
}

*/

func TestReadVcfToAlleleCounts(t *testing.T) {
	samplemap := ReadVcfToAlleleCounts("testdata/human_chrM.vcf")
	if samplemap == nil {
		t.Errorf("Problem with ReadVcfToAlleleCounts")
	}
}

func TestFilterAlleles(t *testing.T) {
	samplemap := ReadVcfToAlleleCounts("testdata/human_chrM.vcf")
	FilterAlleles(samplemap, 100000)
	if len(samplemap) > 0 {
		t.Errorf("Problem with FilterAlleles")
	}
}

func TestFindMajorAllele(t *testing.T) {
	var a int32 = 10
	var c int32 = 50
	var g int32 = 0
	var tt int32 = 2
	var ins int32 = 0
	var del int32 = 1

	major := FindMajorAllele(a, c, g, tt, ins, del)

	if major != 50 {
		t.Errorf("Problem with FindMajorAllele")
	}
}

func TestFindMinorAllele(t *testing.T) {
	var a int32 = 10
	var c int32 = 50
	var g int32 = 0
	var tt int32 = 2
	var ins int32 = 0
	var del int32 = 1

	minor := FindMinorAllele(a, c, g, tt, ins, del)

	if minor != 10 {
		t.Errorf("Problem with FindMinorAllele")
	}
}
