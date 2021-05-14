package vcf

import (
	"strings"
	"testing"
)

var IG1 GenomeSample = GenomeSample{1, 0, false, []string{""}}
var IG2 GenomeSample = GenomeSample{0, 0, false, []string{""}}
var IG3 GenomeSample = GenomeSample{1, 1, false, []string{""}}

var EG1 GenomeSample = GenomeSample{0, 1, false, []string{""}}
var EG2 GenomeSample = GenomeSample{1, 1, false, []string{""}}
var EG3 GenomeSample = GenomeSample{0, 0, false, []string{""}}

var FirstInputSamples []GenomeSample = []GenomeSample{IG1, IG2, IG3}
var FirstExpectedSamples []GenomeSample = []GenomeSample{EG1, EG2, EG3}

var FirstInputVcf Vcf = Vcf{Chr: "chr2", Pos: 4, Ref: "C", Alt: strings.Split("T", ","), Samples: FirstInputSamples}
var FirstExpectedVcf Vcf = Vcf{Chr: "chr2", Pos: 4, Ref: "T", Alt: strings.Split("C", ","), Samples: FirstExpectedSamples}

var InvertVcfTests = []struct {
	Input    Vcf
	Expected Vcf
}{
	{FirstInputVcf, FirstExpectedVcf},
}

func TestInvertVcf(t *testing.T) {
	var answer Vcf
	for _, v := range InvertVcfTests {
		answer = InvertVcf(v.Input)
		if !isEqual(answer, v.Expected) {
			t.Errorf("Error in InvertVCF. Input Ref: %s. Expected Ref: %s.", v.Input.Ref, v.Expected.Ref)
		}
	}
}
