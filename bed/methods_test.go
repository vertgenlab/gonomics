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
}

//note that methods take a special style of function name for testing.
func TestBed_GetChrom(t *testing.T) {
	var currChrom string
	for _, v := range BedMethodTests {
		currChrom = v.b.GetChrom()
		if currChrom != v.expectedChrom {
			t.Errorf("Error in Bed method GetChrom.")
		}
	}
}

func TestBed_GetChromStart(t *testing.T) {
	var currChromStart int
	for _, v := range BedMethodTests {
		currChromStart = v.b.GetChromStart()
		if currChromStart != v.expectedChromStart {
			t.Errorf("Error in Bed method GetChromStart.")
		}
	}
}

func TestBed_GetChromEnd(t *testing.T) {
	var currChromEnd int
	for _, v := range BedMethodTests {
		currChromEnd = v.b.GetChromEnd()
		if currChromEnd != v.expectedChromEnd {
			t.Errorf("Error in Bed method GetChromEnd.")
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
