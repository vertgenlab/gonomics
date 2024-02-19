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

var roundToSigFigsTests = []struct {
	Num		float64
	SigFigs	int
	Expected float32	
}{
	{Num: 0.2837562,
	SigFigs: 4,
	Expected: 0.2838,
		},
	{Num: 0.2837562,
		SigFigs: 7,
		Expected: 0.2837562,
		},
	{Num: 1.2837523,
		SigFigs: 3,
		Expected: 1.28,
		},	
}

func TestRoundToSigFigs(t *testing.T) {
	for _, testCase := range roundToSigFigsTests {
		res := roundToSigFigs(testCase.Num, testCase.SigFigs)
		if res != testCase.Expected {
			t.Errorf("Error: in in pFasta. roundToSigFigs test not as expected.")
		}
	}
}

var getBaseProbsAtPosTests = []struct {
	Base	pDna.Float32Base
	SigFigs	int
	Expected []float32	
}{
	{Base: pDna.Float32Base{
			A: 0.246,
			C: 0.355,
			G: 0.371,
			T: 0.028,
		},
		SigFigs: 2,	
		Expected: []float32{0.25, 0.36, 0.37, 0.028},
	},
	// {Base: pDna.Float32Base{
	// 			A: 0.246,
	// 			C: 0.355,
	// 			G: 0.371,
	// 			T: 0.028,
	// 		},
	// 	SigFigs: 1,
	// 	Expected: []float32{0.2, 0.4, 0.4, 0.03},
	// },
	// {Base: pDna.Float32Base{
	// 			A: 0.24694832,
	// 			C: 0.3,
	// 			G: 0.371,
	// 			T: 0.028388585937,
	// 		},
	// 	SigFigs: 4,
	// 	Expected: []float32{0.2469, 0.3, 0.371, 0.02839},
	// },	
}

func TestGetBaseProbsAtPos(t *testing.T) {
	for _, testCase := range getBaseProbsAtPosTests {
		resA, resC, resG, resT := getBaseProbsAtPos(testCase.Base, testCase.SigFigs)
		
		if (resA != testCase.Expected[0] || resC != testCase.Expected[1] || resG != testCase.Expected[2] || resT != testCase.Expected[3]) {
			fmt.Printf("%v, %v, %v, %v\n", resA, resC, resG, resT)
			
			t.Errorf("Error: in in pFasta. getBaseProbsAtPos test not as expected.")
		}
	}
}
