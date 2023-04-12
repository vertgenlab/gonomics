package vcf

import (
	"strings"
	"testing"

	"github.com/vertgenlab/gonomics/dna"
)

func TestFixVcf(t *testing.T) {

	ref := map[string][]dna.Base{
		"test": {dna.A, dna.C, dna.G, dna.T}}
	vcfTest := Vcf{
		Chr: "test",
		Pos: 2,
		Ref: "C",
		Alt: strings.Split("-", ",")}

	correctAnswer := Vcf{
		Chr: "test",
		Pos: 1,
		Ref: "AC",
		Alt: strings.Split("A", ",")}

	vcfTest = FixVcf(vcfTest, ref)

	if !isEqual(vcfTest, correctAnswer) {
		t.Errorf("ERROR: Problem fixing VCF. Expected %v, got %v", correctAnswer, vcfTest)
	}
}
