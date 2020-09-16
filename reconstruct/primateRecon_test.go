package reconstruct

import (
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/fasta"
	"testing"
)

var hum1 = dna.StringToBases("AATNAAATTTCGTATC")
var bon1 = dna.StringToBases("AATNA--TTTCCTCTG")
var chimp1 = dna.StringToBases("AATNA--TTTCGTCTA")
var gor1 = dna.StringToBases("AATNA--T--CGTCTG")
var oran1 = dna.StringToBases("AATNA--T--CGTCTG")
var Input1 = []*fasta.Fasta{{"hg38", hum1}, {"panPan2", bon1}, {"panTro6", chimp1}, {"gorGor5", gor1}, {"ponAbe3", oran1}}

var QuickPrimateReconTests = []struct {
	records []*fasta.Fasta
	answer  *fasta.Fasta
}{
	{Input1, &fasta.Fasta{"Human_Chimp_Ancestor", dna.StringToBases("AATNA--TNNCGTCTG")}},
}

func TestQuickPrimateRecon(t *testing.T) {
	for _, test := range QuickPrimateReconTests {
		calculated := QuickPrimateRecon(test.records, true)
		if !fasta.IsEqual(test.answer, calculated) {
			t.Errorf("Problem in QuickPrimateRecon. Expected: %s. Calculated: %s.", dna.BasesToString(test.answer.Seq), dna.BasesToString(calculated.Seq))
		}
	}
}
