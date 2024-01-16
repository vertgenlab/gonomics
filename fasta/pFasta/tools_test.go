package pFasta

import (
	"testing"

	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/dna/pDna"
)

// ExtractTestsOne tests a valid input to Extract
var ExtractTestsSuccess = []struct {
	Input			PFasta
	Start			int
	End				int
	Expected	 	PFasta
}{
	{Input: PFasta{Name: "chr1",
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
	Start: 1,
	End: 3,
	Expected:
		PFasta{Name: "chr1",
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
	},
}

// Want to test that function is failing on invalid input
// // ExtractTestsTwo tests an invalid input to Extract: start > end, should fail
// var ExtractTestsTwo = []struct{
// 	Input			PFasta
// 	Start			int
// 	End				int
// }{
// 	{Input: PFasta{
// 		PFasta{Name: "chr1",
// 			Seq: []pDNA.Float32Base{
// 				pDna.Float32Base{
// 					A: 0.2,
// 					C: 0.3,
// 					G: 0.3,
// 					T: 0.2,
// 				},
// 				pDna.Float32Base{
// 					A: 0.1,
// 					C: 0.2,
// 					G: 0.3,
// 					T: 0.4,
// 				},
// 				pDna.Float32Base{
// 					A: 0.25,
// 					C: 0.25,
// 					G: 0.25,
// 					T: 0.25,
// 				},
// 				pDna.Float32Base{
// 					A: 0.6,
// 					C: 0.2,
// 					G: 0.1,
// 					T: 0.1,
// 				},
// 				pDna.Float32Base{
// 					A: 0.2,
// 					C: 0.3,
// 					G: 0.4,
// 					T: 0.1,
// 				},
// 			},
// 		},
// 	},
// 	Start: 4,
// 	End: 3,
// 	},
	
// }

// // ExtractTestsThree tests an invalid input to Extract: start = end, should fail
// var ExtractTestsThree = []struct{
// 	Input			PFasta
// 	Start			int
// 	End				int
// }{
// 	{Input: PFasta{
// 		PFasta{Name: "chr1",
// 			Seq: []pDNA.Float32Base{
// 				pDna.Float32Base{
// 					A: 0.2,
// 					C: 0.3,
// 					G: 0.3,
// 					T: 0.2,
// 				},
// 				pDna.Float32Base{
// 					A: 0.1,
// 					C: 0.2,
// 					G: 0.3,
// 					T: 0.4,
// 				},
// 				pDna.Float32Base{
// 					A: 0.25,
// 					C: 0.25,
// 					G: 0.25,
// 					T: 0.25,
// 				},
// 				pDna.Float32Base{
// 					A: 0.6,
// 					C: 0.2,
// 					G: 0.1,
// 					T: 0.1,
// 				},
// 				pDna.Float32Base{
// 					A: 0.2,
// 					C: 0.3,
// 					G: 0.4,
// 					T: 0.1,
// 					},
// 				},
// 			},
// 		},
// 	},
// 	Start: 3,
// 	End: 3,
// 	Expected: 
// }

// // ExtractTestsFour tests an invalid input to Extract: start out of range, should fail
// var ExtractTestsFour = []struct{
// 	Input			PFasta
// 	Start			int
// 	End				int
// 	Expected	 	PFasta
// }{
// 	{Input: PFasta{
// 		PFasta{Name: "chr1",
// 			Seq: []pDNA.Float32Base{
// 				pDna.Float32Base{
// 					A: 0.2,
// 					C: 0.3,
// 					G: 0.3,
// 					T: 0.2,
// 				},
// 				pDna.Float32Base{
// 					A: 0.1,
// 					C: 0.2,
// 					G: 0.3,
// 					T: 0.4,
// 				},
// 				pDna.Float32Base{
// 					A: 0.25,
// 					C: 0.25,
// 					G: 0.25,
// 					T: 0.25,
// 				},
// 				pDna.Float32Base{
// 					A: 0.6,
// 					C: 0.2,
// 					G: 0.1,
// 					T: 0.1,
// 				},
// 				pDna.Float32Base{
// 					A: 0.2,
// 					C: 0.3,
// 					G: 0.4,
// 					T: 0.1,
// 					},
// 				},
// 			},
// 		},
// 	},
// 	Start: -3,
// 	End: 3,
// 	Expected: 
// }

// // ExtractTestsFive tests an invalid input to Extract: end out of range, should fail
// var ExtractTestsFive = []struct{
// 	Input			PFasta
// 	Start			int
// 	End				int
// 	Expected	 	PFasta
// }{
// 	{Input: PFasta{
// 		PFasta{Name: "chr1",
// 			Seq: []pDNA.Float32Base{
// 				pDna.Float32Base{
// 					A: 0.2,
// 					C: 0.3,
// 					G: 0.3,
// 					T: 0.2,
// 				},
// 				pDna.Float32Base{
// 					A: 0.1,
// 					C: 0.2,
// 					G: 0.3,
// 					T: 0.4,
// 				},
// 				pDna.Float32Base{
// 					A: 0.25,
// 					C: 0.25,
// 					G: 0.25,
// 					T: 0.25,
// 				},
// 				pDna.Float32Base{
// 					A: 0.6,
// 					C: 0.2,
// 					G: 0.1,
// 					T: 0.1,
// 				},
// 				pDna.Float32Base{
// 					A: 0.2,
// 					C: 0.3,
// 					G: 0.4,
// 					T: 0.1,
// 					},
// 				},
// 			},
// 		},
// 	},
// 	Start: 4,
// 	End: 8,
// 	Expected: 
// }

var ExtractBedTestsSuccess = []struct {
	Input			[]PFasta
	Region			bed.Bed
	Expected	 	PFasta
}{
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
	Region: bed.Bed{Chrom: "chr1",
				ChromStart: 3,
				ChromEnd: 	5,
		},
	Expected:
		PFasta{Name: "chr1",
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
					T: 0.1,
				},
			},
		},
	},
}

// var SampleTests = []struct {
// 	Input			PFasta
// }{
// 	{Input: []PFasta{PFasta{Name: "chr1",
// 			Seq: []pDna.Float32Base{
// 				pDna.Float32Base{
// 					A: 0.2,
// 					C: 0.3,
// 					G: 0.3,
// 					T: 0.2,
// 				},
// 				pDna.Float32Base{
// 					A: 0.1,
// 					C: 0.2,
// 					G: 0.3,
// 					T: 0.4,
// 				},
// 				pDna.Float32Base{
// 					A: 0.25,
// 					C: 0.25,
// 					G: 0.25,
// 					T: 0.25,
// 				},
// 				pDna.Float32Base{
// 					A: 0.6,
// 					C: 0.2,
// 					G: 0.1,
// 					T: 0.1,
// 				},
// 				pDna.Float32Base{
// 					A: 0.2,
// 					C: 0.3,
// 					G: 0.4,
// 					T: 0.1,
// 					},
// 				},
// 			},
// 		},
// 	},
// }

func TestExtract(t *testing.T) {
	for _, testCase := range ExtractTestsSuccess {
		res := []PFasta{Extract(testCase.Input, testCase.Start, testCase.End)}
		expect := []PFasta{testCase.Expected}
		if !AllAreEqual(res, expect, 1e-3){
			t.Errorf("Error: in pFasta. Extract valid input test was not as expected.")
		}
	}
}

func TestExtractBed(t *testing.T) {
	for _, testCase := range ExtractBedTestsSuccess {
		res := ExtractBed(testCase.Input, testCase.Region)
		if !Equal(res, testCase.Expected, 1e-3) {
			t.Errorf("Error: in pFasta. ExtractBed valid input test not as expected.\n")
		}
	}
}

// func TestSample(t *testing.T) {
// 	for _, testCase := range SampleTests {
// 		res := make([]fasta.Fasta)
// 		numSamples := 100
// 		for i := 0; i < numSamples; i++ {
// 			res = res.append(res, Sample(testCase.Input)) // res should be a list of numSamples fastas
// 		}
// 		sampledDistrib := PFasta{Name: testCase.Input.Name, Seq: make([]pDNA.Float32Base, len(testCase.Input.seq))}
// 		for _, sample := range res { // iterate through samples to find distribution
// 			for seqIdx, base := range sample { // iterate through 1 sample to check base at each seqIdx, update sampledDistrib at seqIdx
// 				stringBase = dna.BaseToString(base)
// 				if stringBase == "A" {
// 					sampledDistrib.Seq[seqIdx].A++
// 				} else if stringBase == "C" {
// 					sampledDistrib.Seq[seqIdx].C++
// 				} else if stringBase == "G" {
// 					sampledDistrib.Seq[seqIdx].G++
// 				} else {
// 					sampledDistrib.Seq[seqIdx].T++
// 				}
// 			}
// 		}
// 		expected := testCase.Input
// 		for seqIdx, distrib := range sampledDistrib {
// 			sampledDistrib.Seq[seqIdx].A = float32(sampledDistrib[seqIdx].A )/float32(numSamples)
// 			sampledDistrib.Seq[seqIdx].G = float32(sampledDistrib[seqIdx].A )/float32(numSamples)
// 			sampledDistrib.Seq[seqIdx].C = float32(sampledDistrib[seqIdx].A )/float32(numSamples)
// 			sampledDistrib.Seq[seqIdx].T = float32(sampledDistrib[seqIdx].A )/float32(numSamples)
// 		}

// 		if !Equal(sampledDistrib, expected, 1e-3) {
// 			t.Errorf("Error: in pFasta. Sample valid input test not as expected.\n")
// 		}
// 	}
// }
