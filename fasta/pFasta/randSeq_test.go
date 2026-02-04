package pFasta

import (
	"github.com/vertgenlab/gonomics/fileio"
	"math/rand"
	"testing"
)

var RandSeqTests = []struct {
	Length     int
	Name       string
	SeedSet    bool
	SetSeed    int64
	RandSource *rand.Rand
	Expected   string
	OutFile    string
	Precision  float32
}{
	{Length: 5,
		Name:       "chr1",
		SeedSet:    false,
		SetSeed:    7,
		RandSource: nil,
		Expected:   "testdata/randSeq_Expected_0.pfa",
		OutFile:    "testdata/out.randseq.test.pfa",
		Precision:  1e-3,
	},
	{Length: 100,
		Name:       "chr1",
		SeedSet:    true,
		SetSeed:    8,
		RandSource: rand.New(rand.NewSource(2)),
		Expected:   "testdata/randSeq_Expected_1.pfa",
		OutFile:    "testdata/out.randseq.test.pfa",
		Precision:  1e-3,
	},
}

func TestRandSeq(t *testing.T) {
	for _, v := range RandSeqTests {
		observed := []PFasta{RandSeq(v.Length, v.Name, v.SeedSet, v.SetSeed, v.RandSource)}
		Write(v.OutFile, observed)
		expected := Read(v.Expected)
		if !AllAreEqual(expected, observed, v.Precision) {
			t.Errorf("Error: in pFasta. RandSeq test case not as expected.\n")
		} else {
			fileio.MustRemove(v.OutFile)
		}
	}
}
