package pFasta

import (
	"testing"
	"fmt"
	"math/rand"

	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/dna/pDna"
	"github.com/vertgenlab/gonomics/fasta"
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
	{Input: []PFasta{PFasta{Name: "chr3",
		Seq: []pDna.Float32Base{
			pDna.Float32Base{
				A: 0.2,
				C: 0.3,
				G: 0.4,
				T: 0.1,
			},
			pDna.Float32Base{
				A: 0.25,
				C: 0.25,
				G: 0.25,
				T: 0.25,
			},
			pDna.Float32Base{
				A: 0.2,
				C: 0.3,
				G: 0.3,
				T: 0.2,
			},
			pDna.Float32Base{
				A: 0.1,
				C: 0.2,
				G: 0.3,
				T: 0.4,
			},
			pDna.Float32Base{
				A: 0.6,
				C: 0.2,
				G: 0.1,
				T: 0.1,
			},
		},
	},
		PFasta{Name: "chr1",
			Seq: []pDna.Float32Base{
				pDna.Float32Base{
					A: 0.2,
					C: 0.3,
					G: 0.3,
					T: 0.2,
				},
				pDna.Float32Base{
					A: 0.1,
					C: 0.2,
					G: 0.3,
					T: 0.4,
				},
				pDna.Float32Base{
					A: 0.25,
					C: 0.25,
					G: 0.25,
					T: 0.25,
				},
				pDna.Float32Base{
					A: 0.6,
					C: 0.2,
					G: 0.1,
					T: 0.1,
				},
				pDna.Float32Base{
					A: 0.2,
					C: 0.3,
					G: 0.4,
					T: 0.1,
				},
			},
		},
		PFasta{Name: "chr2",
			Seq: []pDna.Float32Base{
				pDna.Float32Base{
					A: 0.2,
					C: 0.3,
					G: 0.4,
					T: 0.1,
				},
				pDna.Float32Base{
					A: 0.25,
					C: 0.25,
					G: 0.25,
					T: 0.25,
				},
				pDna.Float32Base{
					A: 0.2,
					C: 0.3,
					G: 0.3,
					T: 0.2,
				},
				pDna.Float32Base{
					A: 0.1,
					C: 0.2,
					G: 0.3,
					T: 0.4,
				},
				pDna.Float32Base{
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
				pDna.Float32Base{
					A: 0.1,
					C: 0.2,
					G: 0.3,
					T: 0.4,
				},
				pDna.Float32Base{
					A: 0.25,
					C: 0.25,
					G: 0.25,
					T: 0.25,
				},
			},
		},
		Precision: 1e-3,
	},
	{Input: []PFasta{PFasta{Name: "chr1",
		Seq: []pDna.Float32Base{
			pDna.Float32Base{
				A: 0.2,
				C: 0.3,
				G: 0.3,
				T: 0.2,
			},
			pDna.Float32Base{
				A: 0.1,
				C: 0.2,
				G: 0.3,
				T: 0.4,
			},
			pDna.Float32Base{
				A: 0.25,
				C: 0.25,
				G: 0.25,
				T: 0.25,
			},
			pDna.Float32Base{
				A: 0.6,
				C: 0.2,
				G: 0.1,
				T: 0.1,
			},
			pDna.Float32Base{
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
				pDna.Float32Base{
					A: 0.1,
					C: 0.2,
					G: 0.3,
					T: 0.4,
				},
				pDna.Float32Base{
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
	{Input: []PFasta{PFasta{Name: "chr3",
						Seq: []pDna.Float32Base{
							pDna.Float32Base{
								A: 0.2,
								C: 0.3,
								G: 0.4,
								T: 0.1,
							},
							pDna.Float32Base{
								A: 0.25,
								C: 0.25,
								G: 0.25,
								T: 0.25,
							},
							pDna.Float32Base{
								A: 0.2,
								C: 0.3,
								G: 0.3,
								T: 0.2,
							},
							pDna.Float32Base{
								A: 0.1,
								C: 0.2,
								G: 0.3,
								T: 0.4,
							},
							pDna.Float32Base{
								A: 0.6,
								C: 0.2,
								G: 0.1,
								T: 0.1,
							}}},
					PFasta{Name: "chr1",
						Seq: []pDna.Float32Base{
							pDna.Float32Base{
								A: 0.2,
								C: 0.3,
								G: 0.3,
								T: 0.2,
							},
							pDna.Float32Base{
								A: 0.1,
								C: 0.2,
								G: 0.3,
								T: 0.4,
							},
							pDna.Float32Base{
								A: 0.25,
								C: 0.25,
								G: 0.25,
								T: 0.25,
							},
							pDna.Float32Base{
								A: 0.6,
								C: 0.2,
								G: 0.1,
								T: 0.1,
							},
							pDna.Float32Base{
								A: 0.2,
								C: 0.3,
								G: 0.4,
								T: 0.1,
							}}},
					PFasta{Name: "chr2",
						Seq: []pDna.Float32Base{
							pDna.Float32Base{
								A: 0.2,
								C: 0.3,
								G: 0.4,
								T: 0.1,
							},
							pDna.Float32Base{
								A: 0.25,
								C: 0.25,
								G: 0.25,
								T: 0.25,
							},
							pDna.Float32Base{
								A: 0.2,
								C: 0.3,
								G: 0.3,
								T: 0.2,
							},
							pDna.Float32Base{
								A: 0.1,
								C: 0.2,
								G: 0.3,
								T: 0.4,
							},
							pDna.Float32Base{
								A: 0.6,
								C: 0.2,
								G: 0.1,
								T: 0.1,
							}}}},
	Region: []bed.Bed{bed.Bed{Chrom: "chr1", ChromStart: 3, ChromEnd: 5, FieldsInitialized: 3},
				      bed.Bed{Chrom: "chr3", ChromStart: 0, ChromEnd: 3, FieldsInitialized: 3}},
	TakeCoords: false,
	Expected: []PFasta{PFasta{Name: "chr1",
						Seq: []pDna.Float32Base{
							pDna.Float32Base{
								A: 0.6,
								C: 0.2,
								G: 0.1,
								T: 0.1,
							},
							pDna.Float32Base{
								A: 0.2,
								C: 0.3,
								G: 0.4,
								T: 0.1,}}},
					PFasta{Name: "chr3",
						Seq: []pDna.Float32Base{
							pDna.Float32Base{
								A: 0.2,
								C: 0.3,
								G: 0.4,
								T: 0.1,
							},
							pDna.Float32Base{
								A: 0.25,
								C: 0.25,
								G: 0.25,
								T: 0.25,
							},
							pDna.Float32Base{
								A: 0.2,
								C: 0.3,
								G: 0.3,
								T: 0.2,
							}}}},
	Precision: 1e-3},
	{Input: []PFasta{PFasta{Name: "chr3",
						Seq: []pDna.Float32Base{
							pDna.Float32Base{
								A: 0.2,
								C: 0.3,
								G: 0.4,
								T: 0.1,
							},
							pDna.Float32Base{
								A: 0.25,
								C: 0.25,
								G: 0.25,
								T: 0.25,
							},
							pDna.Float32Base{
								A: 0.2,
								C: 0.3,
								G: 0.3,
								T: 0.2,
							},
							pDna.Float32Base{
								A: 0.1,
								C: 0.2,
								G: 0.3,
								T: 0.4,
							},
							pDna.Float32Base{
								A: 0.6,
								C: 0.2,
								G: 0.1,
								T: 0.1,
							}}},
					PFasta{Name: "chr1",
						Seq: []pDna.Float32Base{
							pDna.Float32Base{
								A: 0.2,
								C: 0.3,
								G: 0.3,
								T: 0.2,
							},
							pDna.Float32Base{
								A: 0.1,
								C: 0.2,
								G: 0.3,
								T: 0.4,
							},
							pDna.Float32Base{
								A: 0.25,
								C: 0.25,
								G: 0.25,
								T: 0.25,
							},
							pDna.Float32Base{
								A: 0.6,
								C: 0.2,
								G: 0.1,
								T: 0.1,
							},
							pDna.Float32Base{
								A: 0.2,
								C: 0.3,
								G: 0.4,
								T: 0.1,
							}}},
					PFasta{Name: "chr2",
						Seq: []pDna.Float32Base{
							pDna.Float32Base{
								A: 0.2,
								C: 0.3,
								G: 0.4,
								T: 0.1,
							},
							pDna.Float32Base{
								A: 0.25,
								C: 0.25,
								G: 0.25,
								T: 0.25,
							},
							pDna.Float32Base{
								A: 0.2,
								C: 0.3,
								G: 0.3,
								T: 0.2,
							},
							pDna.Float32Base{
								A: 0.1,
								C: 0.2,
								G: 0.3,
								T: 0.4,
							},
							pDna.Float32Base{
								A: 0.6,
								C: 0.2,
								G: 0.1,
								T: 0.1,
							}}}},
	Region: []bed.Bed{bed.Bed{Chrom: "chr3", ChromStart: 0, ChromEnd: 2, FieldsInitialized: 3},
				      bed.Bed{Chrom: "chr3", ChromStart: 1, ChromEnd: 3, FieldsInitialized: 3}},
	TakeCoords: true,
	Expected: []PFasta{PFasta{Name: "chr3:0-2",
						Seq: []pDna.Float32Base{
							pDna.Float32Base{
								A: 0.2,
								C: 0.3,
								G: 0.4,
								T: 0.1,
							},
							pDna.Float32Base{
								A: 0.25,
								C: 0.25,
								G: 0.25,
								T: 0.25}}},
					PFasta{Name: "chr3:1-3",
						Seq: []pDna.Float32Base{
							pDna.Float32Base{
								A: 0.25,
								C: 0.25,
								G: 0.25,
								T: 0.25,
							},
							pDna.Float32Base{
								A: 0.2,
								C: 0.3,
								G: 0.3,
								T: 0.2}}}},
	Precision: 1e-3},
}

var SampleTests = []struct {
	Input    []PFasta
	Chrom    string
	SetSeed  int64
	Expected fasta.Fasta
}{
	{Input: []PFasta{PFasta{Name: "chr2",
		Seq: []pDna.Float32Base{
			pDna.Float32Base{
				A: 0.2,
				C: 0.3,
				G: 0.4,
				T: 0.1,
			},
			pDna.Float32Base{
				A: 0.25,
				C: 0.25,
				G: 0.25,
				T: 0.25,
			},
			pDna.Float32Base{
				A: 0.2,
				C: 0.3,
				G: 0.3,
				T: 0.2,
			},
			pDna.Float32Base{
				A: 0.1,
				C: 0.2,
				G: 0.3,
				T: 0.4,
			},
			pDna.Float32Base{
				A: 0.6,
				C: 0.2,
				G: 0.1,
				T: 0.1,
			},
		},
	},
		PFasta{Name: "chr1",
			Seq: []pDna.Float32Base{
				pDna.Float32Base{
					A: 0.2,
					C: 0.3,
					G: 0.3,
					T: 0.2,
				},
				pDna.Float32Base{
					A: 0.1,
					C: 0.2,
					G: 0.3,
					T: 0.4,
				},
				pDna.Float32Base{
					A: 0.25,
					C: 0.25,
					G: 0.25,
					T: 0.25,
				},
				pDna.Float32Base{
					A: 0.6,
					C: 0.2,
					G: 0.1,
					T: 0.1,
				},
				pDna.Float32Base{
					A: 0.2,
					C: 0.3,
					G: 0.4,
					T: 0.1,
				},

				pDna.Float32Base{
					A: 0.2,
					C: 0.4,
					G: 0.1,
					T: 0.3,
				},
				pDna.Float32Base{
					A: 0.2,
					C: 0.2,
					G: 0.2,
					T: 0.4,
				},
				pDna.Float32Base{
					A: 0.1,
					C: 0.3,
					G: 0.3,
					T: 0.3,
				},
				pDna.Float32Base{
					A: 0.5,
					C: 0.4,
					G: 0.1,
					T: 0,
				},
				pDna.Float32Base{
					A: 0.1,
					C: 0.1,
					G: 0.7,
					T: 0.1,
				},
			},
		},
		PFasta{Name: "chr3",
			Seq: []pDna.Float32Base{
				pDna.Float32Base{
					A: 0.2,
					C: 0.3,
					G: 0.4,
					T: 0.1,
				},
				pDna.Float32Base{
					A: 0.25,
					C: 0.25,
					G: 0.25,
					T: 0.25,
				},
				pDna.Float32Base{
					A: 0.2,
					C: 0.3,
					G: 0.3,
					T: 0.2,
				},
				pDna.Float32Base{
					A: 0.1,
					C: 0.2,
					G: 0.3,
					T: 0.4,
				},
				pDna.Float32Base{
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

func TestExtract(t *testing.T) {
	for _, testCase := range ExtractTests {
		InFile := "testdata_tools/test_extract_input.pfa"
		Write(InFile, testCase.Input)

		res := Extract(testCase.Input, testCase.Start, testCase.End, testCase.OutputName, testCase.Chrom, testCase.TakeCoords)

		if !IsEqual(res, testCase.Expected, testCase.Precision) {
			t.Errorf("Error: in pFasta. Extract valid input test was not as expected.")
		}

		OutFile := "testdata_tools/test_extract_expected.pfa"
		records := []PFasta{res}
		Write(OutFile, records)
	}
}

func TestExtractBed(t *testing.T) {
	for caseIdx, testCase := range ExtractBedTests {
		InFile := fmt.Sprintf("testdata_tools/test_extractbed_input_%v.pfa", caseIdx)
		Write(InFile, testCase.Input)
		InRegion := fmt.Sprintf("testdata_tools/test_extractbed_input_region_%v.bed", caseIdx)

		bedInput := testCase.Region
		bed.Write(InRegion, bedInput)

		res := ExtractBed(testCase.Input, testCase.Region, testCase.TakeCoords)

		if !AllAreEqual(res, testCase.Expected, testCase.Precision) {
			t.Errorf("Error: in pFasta. ExtractBed valid input test not as expected.\n")
		}

		OutFile := fmt.Sprintf("testdata_tools/test_extractbed_expected_%v.pfa", caseIdx)
		Write(OutFile, res)
	}
}

func TestSample(t *testing.T) {
	for _, testCase := range SampleTests {
		rand.Seed(testCase.SetSeed)
		observed := Sample(testCase.Input, testCase.Chrom)
		if !fasta.IsEqual(observed, testCase.Expected) {
			t.Errorf("Error: in pFasta. Sample valid input test not as expected.\n")
		}
	}
}
