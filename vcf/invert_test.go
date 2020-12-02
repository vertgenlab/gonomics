package vcf

import (
	"testing"
	"strings"
)

var IG1 GenomeSample = GenomeSample{1, 0, false}
var IG2 GenomeSample = GenomeSample{0, 0, false}
var IG3 GenomeSample = GenomeSample{1, 1, false}

var EG1 GenomeSample = GenomeSample{0, 1, false}
var EG2 GenomeSample = GenomeSample{1, 1, false}
var EG3 GenomeSample = GenomeSample{0, 0, false}

var FirstInputGenotypes []GenomeSample = []GenomeSample{IG1, IG2, IG3}
var FirstExpectedGenotypes []GenomeSample = []GenomeSample{EG1, EG2, EG3}

var FirstInputVcf Vcf = Vcf{Chr: "chr2", Pos: 4, Ref: "C", Alt: strings.Split("T", ","), Genotypes: FirstInputGenotypes}
var FirstExpectedVcf Vcf = Vcf{Chr: "chr2", Pos: 4, Ref: "T", Alt: strings.Split("C", ","), Genotypes: FirstExpectedGenotypes}

var InvertVcfTests = []struct {
	Input    Vcf
	Expected Vcf
}{
	{FirstInputVcf, FirstExpectedVcf},
}

func TestInvertVcf(t *testing.T) {
	for _, v := range InvertVcfTests {
		InvertVcf(&v.Input)
		if !isEqual(&v.Input, &v.Expected) {
			t.Errorf("Error in InvertVCF. Input Ref: %s. Expected Ref: %s.", v.Input.Ref, v.Expected.Ref)
		}
	}
}
