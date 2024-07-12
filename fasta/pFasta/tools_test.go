package pFasta

import (
	"math/rand"
	"testing"

	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/dna/pDna"
	"github.com/vertgenlab/gonomics/fasta"
)

var SampleTests = []struct {
	Input    []PFasta
	Chrom    string
	SetSeed  int64
	Expected fasta.Fasta
}{
	{Input: []PFasta{{Name: "chr2",
		Seq: []pDna.Float32Base{
			{
				A: 0.2,
				C: 0.3,
				G: 0.4,
				T: 0.1,
			},
			{
				A: 0.25,
				C: 0.25,
				G: 0.25,
				T: 0.25,
			},
			{
				A: 0.2,
				C: 0.3,
				G: 0.3,
				T: 0.2,
			},
			{
				A: 0.1,
				C: 0.2,
				G: 0.3,
				T: 0.4,
			},
			{
				A: 0.6,
				C: 0.2,
				G: 0.1,
				T: 0.1,
			},
		},
	},
		{Name: "chr1",
			Seq: []pDna.Float32Base{
				{
					A: 0.2,
					C: 0.3,
					G: 0.3,
					T: 0.2,
				},
				{
					A: 0.1,
					C: 0.2,
					G: 0.3,
					T: 0.4,
				},
				{
					A: 0.25,
					C: 0.25,
					G: 0.25,
					T: 0.25,
				},
				{
					A: 0.6,
					C: 0.2,
					G: 0.1,
					T: 0.1,
				},
				{
					A: 0.2,
					C: 0.3,
					G: 0.4,
					T: 0.1,
				},

				{
					A: 0.2,
					C: 0.4,
					G: 0.1,
					T: 0.3,
				},
				{
					A: 0.2,
					C: 0.2,
					G: 0.2,
					T: 0.4,
				},
				{
					A: 0.1,
					C: 0.3,
					G: 0.3,
					T: 0.3,
				},
				{
					A: 0.5,
					C: 0.4,
					G: 0.1,
					T: 0,
				},
				{
					A: 0.1,
					C: 0.1,
					G: 0.7,
					T: 0.1,
				},
			},
		},
		{Name: "chr3",
			Seq: []pDna.Float32Base{
				{
					A: 0.2,
					C: 0.3,
					G: 0.4,
					T: 0.1,
				},
				{
					A: 0.25,
					C: 0.25,
					G: 0.25,
					T: 0.25,
				},
				{
					A: 0.2,
					C: 0.3,
					G: 0.3,
					T: 0.2,
				},
				{
					A: 0.1,
					C: 0.2,
					G: 0.3,
					T: 0.4,
				},
				{
					A: 0.6,
					C: 0.2,
					G: 0.1,
					T: 0.1,
				},
			},
		},
	},
		Chrom:   "chr1",
		SetSeed: 7,
		Expected: fasta.Fasta{Name: "chr1",
			Seq: dna.StringToBases("TCATGACCAG")},
	},
}

func TestSample(t *testing.T) {
	for _, testCase := range SampleTests {
		seed := rand.New(rand.NewSource(testCase.SetSeed))
		observed := Sample(testCase.Input, testCase.Chrom, seed)
		if !fasta.IsEqual(observed, testCase.Expected) {
			t.Errorf("Error: in pFasta. Sample valid input test not as expected.\n")
		}
	}
}
