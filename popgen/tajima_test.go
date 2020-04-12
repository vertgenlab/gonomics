package popgen

import (
	"fmt"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/fasta"
	"testing"
)

var seqA []dna.Base = dna.StringToBases("ATAATAAAAAAATAATAAAAAAATAAAAAAAATAAAAAAA")
var seqB []dna.Base = dna.StringToBases("AAAAAAAATAAATAATAAAAAAATAAAAAAAAAAAAAAAA")
var seqC []dna.Base = dna.StringToBases("AAAATAAAAATATAATAAAAAAATATAAAAAAAAAAAAAA")
var seqD []dna.Base = dna.StringToBases("AAAAAAAAAAAATAATAAAAAAATAAATAAATAAAAAAAA")
var seqE []dna.Base = dna.StringToBases("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA")
var seqF []dna.Base = dna.StringToBases("AAAAAAAAAATATAATAAAAAAATAAAAAAAAAAAAAAAA")
var seqG []dna.Base = dna.StringToBases("AAAAAAAAAATATAAAAAAAAAATAAAAAAAAAATTAAAA")
var seqH []dna.Base = dna.StringToBases("AAAAAAAAAATATAATAAAAAAATAAAAAAAAAAATAAAA")
var testFa []*fasta.Fasta = []*fasta.Fasta{{"eggplant", seqA}, {"raddish", seqB}, {"rhubarb", seqC}, {"asparagus", seqD}, {"broccoli", seqE}, {"tomato", seqF}, {"celery", seqG}, {"carrot", seqH}}
var expected float64 = -1.296575

func TestTajima(t *testing.T) {
	input := Tajima(testFa)
	//fmt.Printf("%f\n",input)
	if fmt.Sprintf("%f", input) != fmt.Sprintf("%f", expected) {
		t.Errorf("Do not match. Input: %f. Expected: %f.", input, expected)
	}
}
