package pFasta

import (
	"fmt"
	"testing"

	"golang.org/x/exp/rand"

	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/dna/pDna"
	"github.com/vertgenlab/gonomics/fasta"
)

// ExtractTests tests a valid input to Extract
var ExtractTests = []struct {
	Input      PFasta
	Start      int
	End        int
	OutputName string
	Expected   PFasta
	Precision  float32
}{
	{Input: PFasta{Name: "chr1",
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
		},
	},
		Start:      1,
		End:        3,
		OutputName: "testChr1",
		Expected: PFasta{Name: "testChr1",
			Seq: []pDna.Float32Base{
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
			},
		},
		Precision: 1e-3,
	},
}

// ExtractBedTests tests a valid input to ExtractBed
var ExtractBedTests = []struct {
	Input      []PFasta
	Region     bed.Bed
	OutputName string
	Expected   PFasta
	Precision  float32
}{
	{Input: []PFasta{{Name: "chr3",
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
			},
		},
		{Name: "chr2",
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
		Region: bed.Bed{Chrom: "chr1",
			ChromStart: 3,
			ChromEnd:   5,
		},
		OutputName: "testChr1",
		Expected: PFasta{Name: "testChr1",
			Seq: []pDna.Float32Base{
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
			},
		},
		Precision: 1e-3,
	},
}

var SampleTests = []struct {
	Input    PFasta
	SetSeed  uint64
	Expected fasta.Fasta
}{
	{Input: PFasta{Name: "chr1",
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
		SetSeed: 7,
		Expected: fasta.Fasta{Name: "chr1",
			Seq: dna.StringToBases("GGAACTCTGG")},
	},
}

func TestExtract(t *testing.T) {
	for _, testCase := range ExtractTests {
		res := Extract(testCase.Input, testCase.Start, testCase.End, testCase.OutputName)
		if !IsEqual(res, testCase.Expected, testCase.Precision) {
			t.Errorf("Error: in pFasta. Extract valid input test was not as expected.")
		}
	}
}

func TestExtractBed(t *testing.T) {
	for _, testCase := range ExtractBedTests {
		res := ExtractBed(testCase.Input, testCase.Region, testCase.OutputName)
		if !IsEqual(res, testCase.Expected, testCase.Precision) {
			t.Errorf("Error: in pFasta. ExtractBed valid input test not as expected.\n")
		}
	}
}

func TestSample(t *testing.T) {
	for _, testCase := range SampleTests {
		rand.Seed(testCase.SetSeed)

		observed := Sample(testCase.Input)
		if !fasta.IsEqual(observed, testCase.Expected) {
			fmt.Println(dna.BasesToString(observed.Seq))
			t.Errorf("Error: in pFasta. Sample valid input test not as expected.\n")
		}

	}
}
