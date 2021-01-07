package reconstruct

import (
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/fasta"
	"testing"
)

var hum2 = dna.StringToBases("AATNA--TTTCGTCTA")
var bon2 = dna.StringToBases("AATNA--TTTCCTCTG")
var chimp2 = dna.StringToBases("AATNAAATTTCGTATC")
var gor2 = dna.StringToBases("AATNA--T--CGTCTG")
var oran2 = dna.StringToBases("AATNA--T--CGTCTG")
var Input2 = []*fasta.Fasta{{"hg38", hum2}, {"panPan2", bon2}, {"panTro6", chimp2}, {"gorGor5", gor2}, {"ponAbe3", oran2}}

var ChimpReconTests = []struct {
	records []*fasta.Fasta
	answer *fasta.Fasta
}{
	{Input2, &fasta.Fasta{"Human_Chimp_Ancestor", dna.StringToBases("AATNA--TNNCGTCTN")}},
}

func TestChimpAncestorRecon(t *testing.T) {
	for _, test := range ChimpReconTests {
		calculated := ChimpAncestorRecon(test.records, true)
		if !fasta.IsEqual(test.answer, calculated) {
			t.Errorf("Problem in ChimpAncestorRecon. Expected: %s. Calculated: %s.", dna.BasesToString(test.answer.Seq), dna.BasesToString(calculated.Seq))
		}
	}
}