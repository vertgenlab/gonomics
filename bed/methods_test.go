package bed

import (
	"testing"
)

var BedStructTests = []struct {
	bed        Bed
	chrom      string
	chromStart int
	chromEnd   int
}{
	{bed: Bed{Chrom: "chr1", ChromStart: 10, ChromEnd: 20, Name: "", Score: 0, Strand: Positive}, chrom: "chr1", chromStart: 10, chromEnd: 20},
	{bed: Bed{Chrom: "chr10", ChromStart: 11, ChromEnd: 18, Name: "", Score: 0, Strand: Positive}, chrom: "chr10", chromStart: 11, chromEnd: 18},
	{bed: Bed{Chrom: "chr12", ChromStart: 12, ChromEnd: 16, Name: "", Score: 0, Strand: Positive}, chrom: "chr12", chromStart: 12, chromEnd: 16},
}

func TestBed_GetChrom(t *testing.T) {
	for _, v := range BedStructTests {
		if v.chrom != v.bed.GetChrom() {
			t.Errorf("Error in bed methods.go GetChrom(); Expected chrom: %v. Actual chrom: %v.", v.chrom, v.bed.GetChrom())
		}
	}
}

func TestBed_GetChromStart(t *testing.T) {
	for _, v := range BedStructTests {
		if v.chromStart != v.bed.GetChromStart() {
			t.Errorf("Error in bed methods.go GetChromStart(); Expected chromStart: %v. Actual chromStart: %v.", v.chromStart, v.bed.GetChromStart())
		}
	}
}

func TestBed_GetChromEnd(t *testing.T) {
	for _, v := range BedStructTests {
		if v.chromEnd != v.bed.GetChromEnd() {
			t.Errorf("Error in bed methods.go GetChromEnd(); Expected chromEnd: %v. Actual chromEnd: %v.", v.chromEnd, v.bed.GetChromEnd())
		}
	}
}






/*
func TestBedSlice_Len(t *testing.T) {
	bedFile := Read("bedFileTest.bed") // this is a type []Bed and I want []*Bed
	if bedFile.Len() != 3 {
		//error here
	}
}
*/

func TestBedSlice_Len(t *testing.T) {
	var allBed BedSlice
	for _, v := range BedStructTests {
		pointer := &v.bed
		allBed = append(allBed, pointer)
		if allBed.Len() != 3 {
			t.Errorf("Error in bed methods.go Len(); Expected %v Actual %v.", 3, allBed.Len())
		}
	}
}
