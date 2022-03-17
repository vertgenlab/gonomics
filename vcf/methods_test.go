package vcf

import "testing"

var VcfMethodsTests = []struct {
	v             Vcf
	expectedChr   string
	expectedStart int
	expectedEnd   int
	newChr        string
	newStart      int
	newEnd        int
}{
	{v: Vcf{Chr: "chr1", Pos: 2, Ref: "A", Alt: []string{"T"}},
		expectedChr:   "chr1",
		expectedStart: 1,
		expectedEnd:   2,
		newChr:        "chr10",
		newStart:      5,
		newEnd:        6,
	},
}

func TestVcf_GetChrom(t *testing.T) {
	for _, v := range VcfMethodsTests {
		if v.v.GetChrom() != v.expectedChr {
			t.Errorf("Error in Vcf method GetChrom.")
		}
	}
}

func TestVcf_GetChromStart(t *testing.T) {
	for _, v := range VcfMethodsTests {
		if v.v.GetChromStart() != v.expectedStart {
			t.Errorf("Error in Vcf method GetChromStart.")
		}
	}
}

func TestVcf_GetChromEnd(t *testing.T) {
	for _, v := range VcfMethodsTests {
		if v.v.GetChromEnd() != v.expectedEnd {
			t.Errorf("Error in Vcf method GetChromStart.")
		}
	}
}

func TestVcf_UpdateCoord(t *testing.T) {
	for _, v := range VcfMethodsTests {
		v.v = v.v.UpdateCoord(v.newChr, v.newStart, v.newEnd).(Vcf)
		if v.v.GetChrom() != v.newChr {
			t.Errorf("Error in Vcf method UpdateCoord, chrom did not match.")
		}
		if v.v.GetChromStart() != v.newStart {
			t.Errorf("Error in Vcf method UpdateCoord, start did not match.")
		}
		if v.v.GetChromEnd() != v.newEnd {
			t.Errorf("Error in Vcf method UpdateCoord, end did not match.")
		}
	}
}
