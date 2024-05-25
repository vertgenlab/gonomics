package vcf

import (
	"testing"
)

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

func TestIsRefWeakAltStrong(t *testing.T) {
	tests := []struct {
		vcf      Vcf
		expected bool
	}{

		{Vcf{Ref: "T", Alt: []string{"C"}}, true},
		{Vcf{Ref: "T", Alt: []string{"G"}}, true},
		{Vcf{Ref: "A", Alt: []string{"C", "G"}}, false},
	}

	for _, test := range tests {
		result := IsRefWeakAltStrong(test.vcf)
		if result != test.expected {
			t.Errorf("IsRefWeakAltStrong(%v) = %v; want %v", test.vcf, result, test.expected)
		}
	}
}

func TestRefStrongAltWeak(t *testing.T) {
	tests := []struct {
		vcf      Vcf
		expected bool
	}{
		{Vcf{Ref: "C", Alt: []string{"A"}}, true},
		{Vcf{Ref: "C", Alt: []string{"T"}}, true},
		{Vcf{Ref: "G", Alt: []string{"C"}}, false},
		{Vcf{Ref: "C", Alt: []string{"A", "T"}}, false},
	}

	for _, test := range tests {
		result := IsRefStrongAltWeak(test.vcf)
		if result != test.expected {
			t.Errorf("Error: IsRefStrongAltWeak(%v) = %v; want %v", test.vcf, result, test.expected)
		}
	}
}

func TestIsNotRefStrongAltWeak(t *testing.T) {
	tests := []struct {
		vcf      Vcf
		expected bool
	}{
		{Vcf{Ref: "C", Alt: []string{"A"}}, false},
		{Vcf{Ref: "G", Alt: []string{"C"}}, true},
	}

	for _, test := range tests {
		result := IsNotRefStrongAltWeak(test.vcf)
		if result != test.expected {
			t.Errorf("Error: IsNotRefStrongAltWeak(%v) = %v; want %v", test.vcf, result, test.expected)
		}
	}
}
func TestIsNotRefWeakAltStrong(t *testing.T) {
	tests := []struct {
		vcf      Vcf
		expected bool
	}{
		{Vcf{Ref: "C", Alt: []string{"A"}}, true},
		{Vcf{Ref: "A", Alt: []string{"C"}}, false},
	}

	for _, test := range tests {
		result := IsNotRefWeakAltStrong(test.vcf)
		if result != test.expected {
			t.Errorf("Error: IsNotRefWeakAltStrong(%v) = %v; want %v", test.vcf, result, test.expected)
		}
	}
}

func TestIsWeakToStrongOrStrongToWeak(t *testing.T) {
	tests := []struct {
		vcf      Vcf
		expected bool
	}{
		{Vcf{Ref: "G", Alt: []string{"T"}}, true},
		{Vcf{Ref: "G", Alt: []string{"C"}}, false},
	}

	for _, test := range tests {
		result := IsWeakToStrongOrStrongToWeak(test.vcf)
		if result != test.expected {
			t.Errorf("Error: IsWeakToStrongOrStrongToWeak(%v) = %v; want %v", test.vcf, result, test.expected)
		}
	}
}

func TestIsNotWeakToStrongOrStrongToWeak(t *testing.T) {
	tests := []struct {
		vcf      Vcf
		expected bool
	}{
		{Vcf{Ref: "A", Alt: []string{"C"}}, false},
		{Vcf{Ref: "A", Alt: []string{"G"}}, false},
		{Vcf{Ref: "T", Alt: []string{"C"}}, false},
		{Vcf{Ref: "T", Alt: []string{"G"}}, false},
	}

	for _, test := range tests {
		result := IsNotWeakToStrongOrStrongToWeak(test.vcf)
		if result != test.expected {
			t.Errorf("Error: IsNotWeakToStrongOrStrongToWeak(%v) = %v; want %v", test.vcf, result, test.expected)
		}
	}
}
