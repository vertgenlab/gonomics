package sam

import (
	"fmt"
	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/cigar"
	"github.com/vertgenlab/gonomics/dna"
	"testing"
)

func TestSamBedToBases(t *testing.T) {
	var obs []dna.Base
	var err error
	rec := []Sam{
		{RName: "chrA", Pos: 101, Cigar: cigar.FromString("10M"), Seq: dna.StringToBases("AACCTTGGAA")},
		{RName: "chrA", Pos: 101, Cigar: cigar.FromString("1S4M2N5M"), Seq: dna.StringToBases("AACCTTGGAA")},
		{RName: "chrA", Pos: 103, Cigar: cigar.FromString("4M2I4M"), Seq: dna.StringToBases("AACCTTGGAA")},
		{RName: "chrA", Pos: 101, Cigar: cigar.FromString("7M5N3M"), Seq: dna.StringToBases("AACCTTGGAA")},
		{RName: "chrA", Pos: 101, Cigar: cigar.FromString("7M1D3M"), Seq: dna.StringToBases("AACCTTGGAA")},
	}
	b := bed.Bed{Chrom: "chrA", ChromStart: 105, ChromEnd: 110}
	exp := []string{"TGGAA", "TGGA", "CTTGGAA", "TG", "TGGA"}
	for i := range rec {
		obs, err = SamBedToBases(rec[i], b)
		if err != nil {
			fmt.Println(err)
		} else if exp[i] != dna.BasesToString(obs) {
			t.Errorf("FAIL! Expected %s, observed: %s", exp[i], dna.BasesToString(obs))
		}
	}
}
