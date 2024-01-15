package pFasta

import (
	"testing"

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

func TestExtract(t *testing.T) {
	for _, testCase := range ExtractTestsSuccess {
		res := []PFasta{Extract(testCase.Input, testCase.Start, testCase.End)}
		expect := []PFasta{testCase.Expected}
		if !AllAreEqual(res, expect, 1e-3){
			t.Errorf("Error: in pFasta. Extract valid input test was not as expected.")
		}
	}
}

