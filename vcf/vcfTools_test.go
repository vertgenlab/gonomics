package vcf

import "testing"

func TestSnp(t *testing.T) {
	tests := []struct {
		vcf      Vcf
		expected bool
	}{
		{Vcf{Info: "SVTYPE=SNP"}, true},
		{Vcf{Info: "SVTYPE=INS"}, false},
		{Vcf{Info: "SVTYPE=DEL"}, false},
		{Vcf{Info: "SOMEOTHERINFO;SVTYPE=SNP;MOREINFO"}, true},
	}

	for _, test := range tests {
		result := Snp(test.vcf)
		if result != test.expected {
			t.Errorf("Snp(%v) = %v; want %v", test.vcf, result, test.expected)
		}
	}
}

func TestIns(t *testing.T) {
	tests := []struct {
		vcf      Vcf
		expected bool
	}{
		{Vcf{Info: "SVTYPE=INS"}, true},
		{Vcf{Info: "SVTYPE=SNP"}, false},
		{Vcf{Info: "SVTYPE=DEL"}, false},
		{Vcf{Info: "SOMEOTHERINFO;SVTYPE=INS;MOREINFO"}, true},
	}

	for _, test := range tests {
		result := Ins(test.vcf)
		if result != test.expected {
			t.Errorf("Ins(%v) = %v; want %v", test.vcf, result, test.expected)
		}
	}
}

func TestDel(t *testing.T) {
	tests := []struct {
		vcf      Vcf
		expected bool
	}{
		{Vcf{Info: "SVTYPE=DEL"}, true},
		{Vcf{Info: "SVTYPE=SNP"}, false},
		{Vcf{Info: "SVTYPE=INS"}, false},
		{Vcf{Info: "SOMEOTHERINFO;SVTYPE=DEL;MOREINFO"}, true},
	}

	for _, test := range tests {
		result := Del(test.vcf)
		if result != test.expected {
			t.Errorf("Del(%v) = %v; want %v", test.vcf, result, test.expected)
		}
	}
}
