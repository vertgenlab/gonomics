package popgen

import (
	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/fasta"
	"testing"
)

var DunnTests = []struct {
	B                    bed.Bed
	Aln                  []fasta.Fasta
	G                    []*Group
	ExpectedDunnValue    float64
	ExpectedSegSiteCount int
	ExpectedMissing      string
}{
	{B: bed.Bed{Chrom: "chr1", ChromStart: 1, ChromEnd: 5, Name: "FirstTest"},
		Aln:                  myAln,
		G:                    myGroups,
		ExpectedDunnValue:    0.5,
		ExpectedSegSiteCount: 3,
		ExpectedMissing:      "Missing: Apple: Tomato: ",
	},
	{B: bed.Bed{Chrom: "chr1", ChromStart: 10, ChromEnd: 20, Name: "SecondTest"},
		Aln:                  myAln,
		G:                    myGroups,
		ExpectedDunnValue:    2,
		ExpectedSegSiteCount: 5,
		ExpectedMissing:      "Missing: Apple: Tomato: ",
	},
}

var myAln []fasta.Fasta = []fasta.Fasta{
	{Name: "honeycrisp", Seq: dna.StringToBases("AACAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA")},
	{Name: "pinkLady", Seq: dna.StringToBases("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA")},
	{Name: "gala", Seq: dna.StringToBases("AATAAAAAAAAAAAAAGAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA")},
	{Name: "sanMarzano", Seq: dna.StringToBases("AAAGAAAAATTTTTAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA")},
	{Name: "roma", Seq: dna.StringToBases("AAAGAAAAATTTGTAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA")},
	{Name: "cherry", Seq: dna.StringToBases("AAATTAAAATTGTTAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA")},
}

var myGroups []*Group = []*Group{
	{Name: "Apple", Members: []string{"honeycrisp", "pinkLady", "gala"}},
	{Name: "Tomato", Members: []string{"sanMarzano", "roma", "cherry"}},
}

func TestDunn(t *testing.T) {
	var testDunnValue float64
	var testSegSiteCount int
	var testMissing string
	for _, v := range DunnTests {
		testDunnValue, testSegSiteCount, testMissing = Dunn(v.B, v.Aln, v.G)
		if testDunnValue != v.ExpectedDunnValue {
			t.Errorf("Error in Dunn. Output Dunn index value was not as expected.")
		}
		if testSegSiteCount != v.ExpectedSegSiteCount {
			t.Errorf("Error in Dunn. Output segSite count not as expeted.")
		}
		if testMissing != v.ExpectedMissing {
			t.Errorf("Error in Dunn. Output MISSING format string not as expected. Found: %v.", testMissing)
		}
	}
}
