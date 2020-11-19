package vcf

import (
	"github.com/vertgenlab/gonomics/dna"
	"testing"
)

var FirstInputVcf Vcf = Vcf{Chr: "chr2", Pos: 4, Ref: "C", Alt: "T"}
var FirstExpectedVcf Vcf = Vcf{Chr: "chr2", Pos: 4, Ref: "T", Alt: "C"}

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

var FirstInputSeq [][]dna.Base = [][]dna.Base{dna.StringToBases("C"), dna.StringToBases("T")}
var FirstExpectedSeq [][]dna.Base = [][]dna.Base{dna.StringToBases("T"), dna.StringToBases("C")}
var IG1 GenomeSample = GenomeSample{1, 0, false}
var IG2 GenomeSample = GenomeSample{0, 0, false}
var IG3 GenomeSample = GenomeSample{1, 1, false}

var EG1 GenomeSample = GenomeSample{0, 1, false}
var EG2 GenomeSample = GenomeSample{1, 1, false}
var EG3 GenomeSample = GenomeSample{0, 0, false}

var FirstInputGenotypes []GenomeSample = []GenomeSample{IG1, IG2, IG3}
var FirstExpectedGenotypes []GenomeSample = []GenomeSample{EG1, EG2, EG3}

var FirstInputGVcf GVcf = GVcf{Vcf: FirstInputVcf, Seq: FirstInputSeq, Genotypes: FirstInputGenotypes}
var FirstExpectedGVcf GVcf = GVcf{Vcf: FirstExpectedVcf, Seq: FirstExpectedSeq, Genotypes: FirstExpectedGenotypes}

var InvertGVcfTests = []struct {
	Input    GVcf
	Expected GVcf
}{
	{FirstInputGVcf, FirstExpectedGVcf},
}

func TestInvertGVcf(t *testing.T) {
	for _, v := range InvertGVcfTests {
		InvertGVcf(&v.Input)
		if !EqualGVcf(v.Input, v.Expected) {
			t.Errorf("Error in TestInvertGVcf.")
		}
	}
}
