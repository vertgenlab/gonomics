package vcf

import (
	"github.com/vertgenlab/gonomics/dna"
	"log"
	"testing"
)

func TestFixVcf(t *testing.T) {

	ref := map[string][]dna.Base{
		"test": []dna.Base{dna.A, dna.C, dna.G, dna.T}}
	vcfTest := Vcf{
		Chr: "test",
		Pos: 2,
		Ref: "C",
		Alt: "-"}

	correctAnswer := Vcf{
		Chr: "test",
		Pos: 1,
		Ref: "AC",
		Alt: "A"}

	FixVcf(&vcfTest, ref)

	if vcfTest != correctAnswer {
		log.Fatalf("ERROR: Problem fixing VCF. Expected %v, got %v", correctAnswer, vcfTest)
	}
}
