package bed

import "testing"

var BedMethodTests = []struct {
	b                  Bed
	expectedChrom      string
	expectedChromStart int
	expectedChromEnd   int
	newChrom           string
	newStart           int
	newEnd             int
}{
	{b: Bed{Chrom: "chr1", ChromStart: 12, ChromEnd: 15},
		expectedChrom:      "chr1",
		expectedChromStart: 12,
		expectedChromEnd:   15,
		newChrom:           "chrX",
		newStart:           50,
		newEnd:             60,
	},
	{b: Bed{Chrom: "chr12", ChromStart: 22, ChromEnd: 30},
		expectedChrom:      "chr12",
		expectedChromStart: 22,
		expectedChromEnd:   30,
		newChrom:           "chr22",
		newStart:           52,
		newEnd:             62,
	},
}

//note that methods take a special style of function name for testing.
func TestBed_GetChrom(t *testing.T) {
	var currChrom string
	for _, v := range BedMethodTests {
		currChrom = v.b.GetChrom()
		if currChrom != v.expectedChrom {
			t.Errorf("Error in Bed method GetChrom. Expected chrom: %v. Actual chrom: %v.", v.expectedChrom, currChrom)
		}
	}
}

func TestBed_GetChromStart(t *testing.T) {
	var currChromStart int
	for _, v := range BedMethodTests {
		currChromStart = v.b.GetChromStart()
		if currChromStart != v.expectedChromStart {
			t.Errorf("Error in Bed method GetChromStart. Expected: %v, Actual: %v", v.expectedChromStart, currChromStart)
		}
	}
}

func TestBed_GetChromEnd(t *testing.T) {
	var currChromEnd int
	for _, v := range BedMethodTests {
		currChromEnd = v.b.GetChromEnd()
		if currChromEnd != v.expectedChromEnd {
			t.Errorf("Error in Bed method GetChromEnd. Expected: %v, Actual %v", v.expectedChromEnd, currChromEnd)
		}
	}
}

func TestBed_UpdateCoord(t *testing.T) {
	for _, v := range BedMethodTests {
		v.b = v.b.UpdateCoord(v.newChrom, v.newStart, v.newEnd).(Bed)
		if v.b.GetChrom() != v.newChrom {
			t.Errorf("Error in Bed method UpdateCoord, chrom did not match.")
		}
		if v.b.GetChromStart() != v.newStart {
			t.Errorf("Error in Bed method UpdateCoord, start did not match.")
		}
		if v.b.GetChromEnd() != v.newEnd {
			t.Errorf("Error in Bed method UpdateCoord, end did not match.")
		}
	}
}

//	Tests for methods that use: type BedSlice []*Bed
func TestBedSlice_Len(t *testing.T) {
	var allBed BedSlice
	for _, v := range BedMethodTests {
		pointer := &v.b //must create a pointer because type BedSlice is []*Bed
		allBed = append(allBed, pointer)
	}
	if allBed.Len() != 2 {
		t.Errorf("Error in bed methods.go Len(); Expected: %v Actual: %v.", 2, allBed.Len())
	}
}

func TestBedSlice_Swap(t *testing.T) {
	var allBed BedSlice
	expectedChrom := BedMethodTests[1].b.Chrom
	for _, v := range BedMethodTests {
		pointer := &v.b //must create a pointer because type BedSlice is []*Bed
		allBed = append(allBed, pointer)
	}
	allBed.Swap(0, 1)
	if allBed[1].Chrom != expectedChrom {
		t.Errorf("Error in bed methods.go Swap(); Expected: %v, Actual: %v.", expectedChrom, allBed[1].Chrom)
	}
}
