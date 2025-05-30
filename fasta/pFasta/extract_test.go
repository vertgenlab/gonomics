package pFasta

import (
	"fmt"
	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/dna/pDna"
	"testing"
)

// ExtractTests tests a valid input to Extract
var ExtractTests = []struct {
	Input      []PFasta
	Start      int
	End        int
	OutputName string
	Chrom      string
	TakeCoords bool
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
		Start:      1,
		End:        3,
		OutputName: "testChr1",
		Chrom:      "chr1",
		TakeCoords: false,
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
	{Input: []PFasta{{Name: "chr1",
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
	},
		Start:      1,
		End:        3,
		OutputName: "testChr1",
		Chrom:      "chr1",
		TakeCoords: false,
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
	Region     []bed.Bed
	TakeCoords bool
	Expected   []PFasta
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
			}}},
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
				}}},
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
				}}}},
		Region: []bed.Bed{{Chrom: "chr1", ChromStart: 3, ChromEnd: 5, FieldsInitialized: 3},
			{Chrom: "chr3", ChromStart: 0, ChromEnd: 3, FieldsInitialized: 3}},
		TakeCoords: false,
		Expected: []PFasta{{Name: "chr1",
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
					T: 0.1}}},
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
					}}}},
		Precision: 1e-3},
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
			}}},
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
				}}},
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
				}}}},
		Region: []bed.Bed{{Chrom: "chr3", ChromStart: 0, ChromEnd: 2, FieldsInitialized: 3},
			{Chrom: "chr3", ChromStart: 1, ChromEnd: 3, FieldsInitialized: 3}},
		TakeCoords: true,
		Expected: []PFasta{{Name: "chr3:0-2",
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
					T: 0.25}}},
			{Name: "chr3:1-3",
				Seq: []pDna.Float32Base{
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
						T: 0.2}}}},
		Precision: 1e-3},
}

func TestExtract(t *testing.T) {
	for _, testCase := range ExtractTests {
		InFile := "testdata/test_extract_input.pfa"
		Write(InFile, testCase.Input)

		res := Extract(testCase.Input, testCase.Start, testCase.End, testCase.OutputName, testCase.Chrom, testCase.TakeCoords)

		if !IsEqual(res, testCase.Expected, testCase.Precision) {
			t.Errorf("Error: in pFasta. Extract valid input test was not as expected.")
		}

		OutFile := "testdata/test_extract_expected.pfa"
		records := []PFasta{res}
		Write(OutFile, records)
	}
}

func TestExtractBed(t *testing.T) {
	for caseIdx, testCase := range ExtractBedTests {
		InFile := fmt.Sprintf("testdata/test_extractbed_input_%v.pfa", caseIdx)
		Write(InFile, testCase.Input)
		InRegion := fmt.Sprintf("testdata/test_extractbed_input_region_%v.bed", caseIdx)

		bedInput := testCase.Region
		bed.Write(InRegion, bedInput)

		res := ExtractBed(testCase.Input, testCase.Region, testCase.TakeCoords)

		if !AllAreEqual(res, testCase.Expected, testCase.Precision) {
			t.Errorf("Error: in pFasta. ExtractBed valid input test not as expected.\n")
		}

		OutFile := fmt.Sprintf("testdata/test_extractbed_expected_%v.pfa", caseIdx)
		Write(OutFile, res)
	}
}
