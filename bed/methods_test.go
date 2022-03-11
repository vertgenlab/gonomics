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

func TestBed_UpdateLift(t *testing.T) {
	for _, v := range BedStructTests {
		v.bed.UpdateLift("chr4", 3, 4)
		if (v.bed.Chrom != "chr4") && (v.bed.ChromStart != 3) && (v.bed.ChromEnd != 4) {
			t.Errorf("Error in bed methods.go UpdateLift(), Expected chr4, 3, 4 Actual: %v, %v, %v", v.bed.Chrom, v.bed.ChromStart, v.bed.ChromEnd)
		}
	}
}

func TestBedSlice_Len(t *testing.T) {
	var allBed BedSlice
	for _, v := range BedStructTests {
		pointer := &v.bed
		allBed = append(allBed, pointer)
	}
	if allBed.Len() != 3 {
		t.Errorf("Error in bed methods.go Len(); Expected: %v Actual: %v.", 3, allBed.Len())
	}
}

func TestBedSlice_Swap(t *testing.T) {
	var allBed BedSlice
	for _, v := range BedStructTests {
		pointer := &v.bed
		allBed = append(allBed, pointer)
	}
	allBed.Swap(1, 2)
	if allBed[1].Chrom != "chr12" {
		t.Errorf("Error in bed methods.go Swap(); Expected: chr12 Actual: %v.", allBed[1].Chrom)
	}
}
