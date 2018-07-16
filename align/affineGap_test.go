package align

import (
	"fmt"
	"github.com/craiglowe/gonomics/common"
	"github.com/craiglowe/gonomics/dna"
	"github.com/craiglowe/gonomics/fasta"
	"testing"
)

var affineAlignTests = []struct {
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
}

var affineAlignChunkTests = []struct {
	seqOne string
	seqTwo string
	aln    string
}{
	{"ACG", "ACG", "ACG\nACG\n"},
	{"ACG", "CCG", "ACG\nCCG\n"},
	{"TTGTTCTTCTTCTTC", "TTGTTCTTCTTATTATTATTCTTC", "TTGTTCTTC---------TTCTTC\nTTGTTCTTCTTATTATTATTCTTC\n"},
	{"ACAACAATAAGAAAAACAAAA", "ACAACAAAAACAAAA", "ACAACAATAAGAAAAACAAAA\nACAACA------AAAACAAAA\n"},
}

var affineScoreTests = []struct {
	filename string
}{
	{"testdata/tmp.aln.fa"},
	{"testdata/hand.aln.fa"},
}

func TestAffineGap(t *testing.T) {
	for _, test := range affineAlignTests {
		basesOne, _ := dna.StringToBases(test.seqOne)
		basesTwo, _ := dna.StringToBases(test.seqTwo)
		_, cigar := AffineGap(basesOne, basesTwo, DefaultScoreMatrix, -400, -30)
		prettyAlignment := View(basesOne, basesTwo, cigar)
		if prettyAlignment != test.aln {
			t.Errorf("The affine gap alignment of %s and %s gave %s, but %s was expected", test.seqOne, test.seqTwo, prettyAlignment, test.aln)
		}
	}
}

func TestAffineGapChunk(t *testing.T) {
	for _, test := range affineAlignChunkTests {
		basesOne, _ := dna.StringToBases(test.seqOne)
		basesTwo, _ := dna.StringToBases(test.seqTwo)
		_, cigar := AffineGapChunk(basesOne, basesTwo, DefaultScoreMatrix, -400, -30, 3)
		prettyAlignment := View(basesOne, basesTwo, cigar)
		if prettyAlignment != test.aln {
			t.Errorf("The affine gap chunk alignment of %s and %s gave %s, but %s was expected", test.seqOne, test.seqTwo, prettyAlignment, test.aln)
		}
	}
}

func TestAffineGapMulti(t *testing.T) {
	for _, test := range affineAlignTests {
		basesOne, _ := dna.StringToBases(test.seqOne)
		basesTwo, _ := dna.StringToBases(test.seqTwo)
		one := []fasta.Fasta{{Name: "one", Seq: basesOne}}
		two := []fasta.Fasta{{Name: "two", Seq: basesTwo}}
		_, cigar := multipleAffineGap(one, two, DefaultScoreMatrix, -400, -30)
		answer := mergeMultipleAlignments(one, two, cigar)
		pretty := fmt.Sprintf("%s\n%s\n", dna.BasesToString(answer[0].Seq), dna.BasesToString(answer[1].Seq))
		if pretty != test.aln {
			t.Errorf("The affine gap alignment of %s and %s gave %s, but %s was expected.  Cigar was: %v", test.seqOne, test.seqTwo, pretty, test.aln, cigar)
		}
	}
}

func TestAffineScore(t *testing.T) {
	for _, test := range affineScoreTests {
		aln, err := fasta.Read(test.filename)
		common.ExitIfError(err)
		score, err := scoreAffineAln(aln[0], aln[1], DefaultScoreMatrix, -10, -1)
		common.ExitIfError(err)
		fmt.Printf("Score of alignment in %s is: %d\n", test.filename, score)
	}
}
