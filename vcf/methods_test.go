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

func TestVcf_String(t *testing.T) {
	tests := []struct {
		record Vcf
		want   string
	}{
		{
			record: Vcf{
				Chr:    "chr1",
				Pos:    123456,
				Id:     "rs123456",
				Ref:    "A",
				Alt:    []string{"C", "G"},
				Qual:   99.9,
				Filter: "PASS",
				Info:   "DP=100",
			},
			want: "chr1\t123456\trs123456\tA\tC,G\t99.9\tPASS\tDP=100\n",
		},
		{
			record: Vcf{
				Chr:    "chr2",
				Pos:    234567,
				Id:     "rs234567",
				Ref:    "G",
				Alt:    []string{"A"},
				Qual:   50.5,
				Filter: "q10",
				Info:   "AF=0.5",
			},
			want: "chr2\t234567\trs234567\tG\tA\t50.5\tq10\tAF=0.5\n",
		},
	}

	for _, tt := range tests {
		t.Run(tt.record.Id, func(t *testing.T) {
			got := tt.record.String()
			if got != tt.want {
				t.Errorf("Error: String() = %q, want %q", got, tt.want)
			}
		})
	}
}

func TestSample_String(t *testing.T) {
	tests := []struct {
		sample   Sample
		expected string
	}{
		{
			sample:   Sample{Alleles: []int16{0, 1}, Phase: []bool{false, true}, FormatData: []string{"GT", "DP"}},
			expected: "0|1:GT:DP",
		},
		{
			sample:   Sample{Alleles: nil, Phase: nil, FormatData: nil},
			expected: ".",
		},
	}
	for _, s := range tests {
		result := s.sample.String()
		if result != s.expected {
			t.Errorf("Error: String() = %q, want %q", result, s.expected)
		}

	}
}

func TestDoubleDigitAlleles(t *testing.T) {
	tests := []struct {
		sample   Sample
		expected string
	}{
		{
			sample:   Sample{Alleles: []int16{1, 50, 20, 100, 999}, Phase: []bool{false, true, false, false, true}, FormatData: []string{"GT", "DP"}},
			expected: "1|50/20/100|999:GT:DP",
		},
		{
			sample:   Sample{Alleles: []int16{999}, Phase: []bool{false}, FormatData: []string{"GT"}},
			expected: "999:GT",
		},
		{
			sample:   Sample{Alleles: []int16{1, 50, 20, 100, 999, 848, 2343, 453}, Phase: []bool{false, true, false, false, true, false, true, true}, FormatData: []string{"GT", "DP"}},
			expected: "1|50/20/100|999/848|2343|453:GT:DP",
		},
	}
	for _, s := range tests {
		result := s.sample.String()
		if result != s.expected {
			t.Errorf("Error: String() = %q, want %q", result, s.expected)
		}

	}
}
