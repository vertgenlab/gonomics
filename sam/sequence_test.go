package sam

import (
	"fmt"
	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/cigar"
	"github.com/vertgenlab/gonomics/dna"
	"testing"
)

func TestSamBedToBases(t *testing.T) {
	b := bed.Bed{Chrom: "enh4.3e_RevComp_ACCAA", ChromStart: 170, ChromEnd: 175}
	s := Sam{RName: "enh4.3e_RevComp_ACCAA", Cigar: cigar.FromString("1S19M56N55M223N76M"), Pos: 151, Seq: dna.StringToBases("ATGATCTAGAGCATGCACCGGGATAAGTGTCGTTGGGGCAGGTCCTCGGGAGAGTTCAGGTTGGTGGGTCCTGGGGCGGGGCCAGGGCGGGGCGTTGGCTATGTCGTAGCACGTGGCCAGGCGCTGCTCGGACTCTGGGAGGCGGAGCTTA")}
	//s := Sam{RName: "enh4.3e_RevComp_ACCAA", Cigar: cigar.FromString("10M5N10M"), Pos: 161, Seq: dna.StringToBases("AAAAAAAAAAGGGGGGGGGG")}
	exp := "ACCAA"
	obs, err := SamBedToBases(s, b)
	if err != nil {
		fmt.Println(err)
	} else if exp != dna.BasesToString(obs) {
		t.Errorf("FAIL! Expected %s, observed: %s", exp, obs)
	}
}
