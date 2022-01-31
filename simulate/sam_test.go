package simulate

import (
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/sam"
	"math/rand"
	"testing"
)

func TestGenerateSam(t *testing.T) {
	rand.Seed(1)
	ref := fasta.Read("testdata/eng.fa")
	actual := IlluminaSam(ref[0].Name, ref[0].Seq, 100)
	expected, _ := sam.Read("testdata/expected.sam")

	if len(actual) != len(expected) {
		t.Error("problem simulating sam file")
	}

	for i := range actual {
		if !sam.Equal(actual[i], expected[i]) {
			t.Error("problem simulating sam file")
		}
	}
}
