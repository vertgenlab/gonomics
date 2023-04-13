package align

import (
	"testing"

	"github.com/vertgenlab/gonomics/fasta"
)

var multiAlignTests = []struct {
	input    string
	expected string
}{
	{"testdata/multiAlignTest.in.fa", "testdata/multiAlignTest.expected.fa"},
	{"testdata/multiAlignTest.in2.fa", "testdata/multiAlignTest.expected2.fa"},
}

func TestMultiAlignGap(t *testing.T) {
	for _, test := range multiAlignTests {
		input := fasta.Read(test.input)

		expected := fasta.Read(test.expected)

		aligned := AllSeqAffine(input, DefaultScoreMatrix, -400, -30)
		alignedChunk := AllSeqAffineChunk(input, DefaultScoreMatrix, -400, -30, 2)

		//fasta.Write("testdata/multiAlignTest.tmp", aligned)

		if !fasta.AllAreEqualIgnoreOrder(aligned, expected) {
			fasta.Write("testdata/multiAlignTest.tmp", aligned)
			t.Errorf("Alignment not as expected: testdata/multiAlignTest.tmp does not equal %s", test.expected)
		}

		if !fasta.AllAreEqualIgnoreOrder(alignedChunk, expected) {
			fasta.Write("testdata/multiAlignTest.tmp", alignedChunk)
			t.Errorf("Alignment not as expected: testdata/multiAlignTest.tmp does not equal %s", test.expected)
		}

	}
}
