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
			name:     "Null FormatData",
			samples:  []Sample{{Alleles: []int16{0, 1}, Phase: []bool{false, false}, FormatData: nil}},
			expected: ".",
		},
		{
			name:     "Null Alleles",
			samples:  []Sample{{Alleles: nil, Phase: []bool{false, false}, FormatData: nil}},
			expected: ".",
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			result := SamplesToString(tt.samples)
			if result != tt.expected {
				t.Errorf("Error: expected %q, got %q", tt.expected, result)
			}
		})
	}
}

func TestPhasedToByte(t *testing.T) {
	tests := []struct {
		input    bool
		expected byte
	}{
		{true, '|'},
		{false, '/'},
	}

	for _, test := range tests {
		result := PhasedToByte(test.input)
		if result != test.expected {
			t.Errorf("Error: PhasedToByte(%v) = %v; want %v", test.input, result, test.expected)
		}
	}
}

// Test for PrintSampleNames function
func TestPrintSampleNames(t *testing.T) {
	header := Header{
		Text: []string{
			"##fileformat=VCFv4.3",
			"##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">",
			"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSample1\tSample2\tSample3",
		},
	}

	expected := "Sample1\nSample2\nSample3\n"
	result := PrintSampleNames(header)

	if result != expected {
		t.Errorf("expected %q, got %q", expected, result)
	}
}
