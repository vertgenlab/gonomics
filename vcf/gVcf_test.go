package vcf

import (
	"testing"
)

func TestSamplesToString(t *testing.T) {
	tests := []struct {
		name     string
		samples  []Sample
		expected string
	}{
		{
			name: "Single sample",
			samples: []Sample{
				{Alleles: []int16{0, 1}, Phase: []bool{false, false}, FormatData: []string{""}},
			},
			expected: "0/1",
		},
		{
			name: "Multiple samples",
			samples: []Sample{
				{Alleles: []int16{0, 1}, Phase: []bool{false, true}, FormatData: []string{"GT", "DP"}},
				{Alleles: []int16{1, 1}, Phase: []bool{true, true}, FormatData: []string{"GT", "DP"}},
			},
			expected: "0|1GT:DP\t1|1GT:DP",
		},
		{
			name:     "No FormatData",
			samples:  nil,
			expected: "",
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			result := SamplesToString(tt.samples)
			if result != tt.expected {
				t.Errorf("expected %q, got %q", tt.expected, result)
			}
		})
	}
}
