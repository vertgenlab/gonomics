package bed

import (
	"testing"
)

var BedStructTests = []struct {
	A                Bed
	chrom			string
	chromStart		int
	chromEnd		int
}{
	{A: Bed{Chrom: "chr1", ChromStart: 10, ChromEnd: 20, Name: "", Score:0, Strand:Positive}, chrom: "chr1", chromStart: 10, chromEnd: 20},
	{A: Bed{Chrom: "chr10", ChromStart: 11, ChromEnd: 18, Name: "", Score:0, Strand:Positive}, chrom: "chr10", chromStart: 11, chromEnd: 18},
	{A: Bed{Chrom: "chr12", ChromStart: 12, ChromEnd: 16, Name: "", Score:0, Strand:Positive}, chrom: "chr12", chromStart: 11, chromEnd: 16},
}


func TestBed_GetChrom(t *testing.T) {
	for _, v := range BedStructTests {
		if v.chrom != v.A.GetChrom() {
			t.Errorf("Error in bed methods.go GetChrom(); Expected chrom: %v. Actual chrom: %v.", v.chrom, v.A.GetChrom())
		}
	}
}


func TestBed_GetChromStart(t *testing.T) {
	for _, v := range BedStructTests {
		if v.chromStart != v.A.GetChromStart() {
			t.Errorf("Error in bed methods.go GetChromStart(); Expected chromStart: %v. Actual chromStart: %v.", v.chrom, v.A.GetChrom())
		}
	}
}

func TestBed_GetChromEnd(t *testing.T) {

}





func TestBedSlice_Len(t *testing.T) {

}

 */