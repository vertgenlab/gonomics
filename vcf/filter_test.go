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

func TestIsSubstitution(t *testing.T) {
	tests := []struct {
		vcf      Vcf
		expected bool
	}{
		{Vcf{Ref: "A", Alt: []string{"C"}}, true},
		{Vcf{Ref: "AT", Alt: []string{"C"}}, false},
		{Vcf{Ref: "A", Alt: []string{"CG"}}, false},
	}

	for _, test := range tests {
		result := IsSubstitution(test.vcf)
		if result != test.expected {
			t.Errorf("Error: IsSubstitution(%v) = %v; want %v", test.vcf, result, test.expected)
		}
	}
}

func TestIsSegregating(t *testing.T) {
	tests := []struct {
		vcf      Vcf
		expected bool
	}{
		{Vcf{Samples: []Sample{}}, false},
		{Vcf{Samples: []Sample{{Alleles: []int16{0, 0}}}}, false},
		{Vcf{Samples: []Sample{{Alleles: []int16{0, 1}}}}, true},
		{Vcf{Samples: []Sample{{Alleles: []int16{1, 1}}}}, false},
		{Vcf{Samples: []Sample{{Alleles: []int16{1}}, {Alleles: []int16{2}}}}, true},
	}

	for _, test := range tests {
		result := IsSegregating(test.vcf)
		if result != test.expected {
			t.Errorf("Error: IsSegregating(%v) = %v; want %v", test.vcf, result, test.expected)
		}
	}
}
