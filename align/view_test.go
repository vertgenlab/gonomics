package align

import (
	"testing"

	"github.com/vertgenlab/gonomics/dna"
)

var alignTests = []struct {
	seqOne string
	seqTwo string
	aln    string
}{
	{"ACGT", "ACGT", "ACGT\nACGT\n"},
	{"ACGT", "CGT", "ACGT\n-CGT\n"},
	{"ACGT", "ACG", "ACGT\nACG-\n"},
	{"CGT", "ACGT", "-CGT\nACGT\n"},
	{"ACG", "ACGT", "ACG-\nACGT\n"},
	{"AGT", "ACGT", "A-GT\nACGT\n"},
	{"ACT", "ACGT", "AC-T\nACGT\n"},
	{"CGCGCGCGCG", "CGCGCGTTTTCGCG", "CGCGCG----CGCG\nCGCGCGTTTTCGCG\n"},
	{"CGCGCGCGCG", "CGAAAACGCGTTTTCGCG", "CG----CGCG----CGCG\nCGAAAACGCGTTTTCGCG\n"},
	{"AA", "GGGAATT", "---AA--\nGGGAATT\n"},
	{"GGGAATT", "AA", "GGGAATT\n---AA--\n"},
	{"AGTACGT", "ACGTACG", "A-GTACGT\nACGTACG-\n"},
}

func TestConstGap(t *testing.T) {
	for _, test := range alignTests {
		basesOne := dna.StringToBases(test.seqOne)
		basesTwo := dna.StringToBases(test.seqTwo)
		_, cigar := ConstGap(basesOne, basesTwo, DefaultScoreMatrix, -430)
		prettyAlignment := View(basesOne, basesTwo, cigar)
		if prettyAlignment != test.aln {
			t.Errorf("The alignment of %s and %s gave\n%s\n, but this was expected\n%s\n", test.seqOne, test.seqTwo, prettyAlignment, test.aln)
		}
	}
}
