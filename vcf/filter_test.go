package vcf

import "testing"

func TestIsHeterozygous(t *testing.T) {
	tests := []struct {
		sample   Sample
		expected bool
	}{
		{Sample{Alleles: []int16{}}, false},
		{Sample{Alleles: []int16{0, 1}}, true},
		{Sample{Alleles: []int16{1, 1}}, false},
	}

	for _, test := range tests {
		result := IsHeterozygous(test.sample)
		if result != test.expected {
			t.Errorf("Error: IsHeterozygous(%v) = %v; want %v", test.sample, result, test.expected)
		}
	}
}

func TestIsHomozygous(t *testing.T) {
	tests := []struct {
		sample   Sample
		expected bool
	}{
		{Sample{Alleles: []int16{}}, false},
		{Sample{Alleles: []int16{1, 1}}, true},
		{Sample{Alleles: []int16{1, 2}}, false},
	}

	for _, test := range tests {
		result := IsHomozygous(test.sample)
		if result != test.expected {
			t.Errorf("Error: IsHomozygous(%v) = %v; want %v", test.sample, result, test.expected)
		}
	}
}

func TestIsBiallelic(t *testing.T) {
	tests := []struct {
		vcf      Vcf
		expected bool
	}{
		{Vcf{Alt: []string{"A"}}, true},
		{Vcf{Alt: []string{"A", "C"}}, false},
	}

	for _, test := range tests {
		result := IsBiallelic(test.vcf)
		if result != test.expected {
			t.Errorf("Error: IsBiallelic(%v) = %v; want %v", test.vcf, result, test.expected)
		}
	}
}
