package browser

import (
	"fmt"
	"github.com/vertgenlab/gonomics/exception"
	"testing"

	"github.com/vertgenlab/gonomics/dna/pDna"
	"github.com/vertgenlab/gonomics/fasta/pFasta"
	"github.com/vertgenlab/gonomics/fileio"
)

var getBaseProbsAtPosTests = []struct {
	Base     pDna.Float32Base
	SigFigs  int
	Expected []float32
}{
	{Base: pDna.Float32Base{
		A: 0.246,
		C: 0.355,
		G: 0.371,
		T: 0.028,
	},
		SigFigs:  2,
		Expected: []float32{0.25, 0.35, 0.37, 0.028}, // want this to be .36
	},
	{Base: pDna.Float32Base{
		A: 0.246,
		C: 0.355,
		G: 0.371,
		T: 0.028,
	},
		SigFigs:  1,
		Expected: []float32{0.2, 0.4, 0.4, 0.03},
	},
	{Base: pDna.Float32Base{
		A: 0.24694837,
		C: 0.3,
		G: 0.371,
		T: 0.02838856,
	},
		SigFigs:  4,
		Expected: []float32{0.2469, 0.3, 0.371, 0.02839},
	},
}

func TestGetBaseProbsAtPos(t *testing.T) {
	for _, testCase := range getBaseProbsAtPosTests {
		resA, resC, resG, resT := getBaseProbsAtPos(testCase.Base, testCase.SigFigs)

		if resA != testCase.Expected[0] || resC != testCase.Expected[1] || resG != testCase.Expected[2] || resT != testCase.Expected[3] {
			fmt.Printf("%v, %v, %v, %v\n", resA, resC, resG, resT)

			t.Errorf("Error: in in pFasta. getBaseProbsAtPos test not as expected.")
		}
	}
}

var PrintOneSetLinesTests = []struct {
	LineLength    int
	SetOfLinesIdx int
	NumIters      int
	LineA         []float32
	LineC         []float32
	LineG         []float32
	LineT         []float32
	Start         int
	Records       []pFasta.PFasta
	Out           string
	SigFigs       int
	LongestName   int
	Expected      string
}{{LineLength: 10,
	SetOfLinesIdx: 0,
	NumIters:      7,
	LineA:         make([]float32, 30),
	LineC:         make([]float32, 30),
	LineG:         make([]float32, 30),
	LineT:         make([]float32, 30),
	Start:         2,
	Records: []pFasta.PFasta{
		pFasta.PFasta{Name: "chr1",
			Seq: []pDna.Float32Base{
				pDna.Float32Base{
					A: 0.23857,
					C: 0.3323,
					G: 0.44958,
					T: 0.139448,
				},
				pDna.Float32Base{
					A: 0.334,
					C: 0.44239,
					G: 0.134875,
					T: 0.12394,
				},
				pDna.Float32Base{
					A: 0.384398,
					C: 0.59723,
					G: 0.3289,
					T: 0.8325,
				},
				pDna.Float32Base{
					A: 0.488532,
					C: 0.23858,
					G: 0.548523,
					T: 0.5473,
				},
				pDna.Float32Base{
					A: 0.92323,
					C: 0.237,
					G: 0.1747,
					T: 0.5839,
				},
				pDna.Float32Base{
					A: 0.483284,
					C: 0.123,
					G: 0.25388,
					T: 0.98243,
				},
				pDna.Float32Base{
					A: 0.00834,
					C: 0.5288,
					G: 0.58001,
					T: 0.4892,
				},
				pDna.Float32Base{
					A: 0.5688,
					C: 0.278,
					G: 0.1937,
					T: 0.1111,
				},
				pDna.Float32Base{
					A: 0.42397,
					C: 0.2358,
					G: 0.984,
					T: 0.3823,
				},
				pDna.Float32Base{
					A: 0.042387,
					C: 0.587,
					G: 0.32871,
					T: 0.73224,
				}}},
	},
	Out:         "testdata/pfa_printOneSetLines_output_toy_1.txt",
	SigFigs:     3,
	LongestName: 15,
	Expected:    "testdata/pfa_printOneSetLines_expected_toy_1.txt",
},
}

// printOneSetLines(lineLength int, setOfLinesIdx int, numIters int, lineA []float32, lineC []float32, lineG []float32, lineT []float32, start int, records []PFasta, out *fileio.EasyWriter, sigFigs int) {

func TestPrintOneSetLines(t *testing.T) {
	for _, testCase := range PrintOneSetLinesTests {
		pFaOut := "testdata/pfa_printOneSetLines_input_toy_1.pfa"
		pFasta.Write(pFaOut, testCase.Records)

		testOut := fileio.EasyCreate(testCase.Out)

		printOneSetLines(testCase.LineLength, testCase.SetOfLinesIdx, testCase.NumIters,
			testCase.LineA, testCase.LineC, testCase.LineG, testCase.LineT, testCase.Start,
			testCase.Records, testOut, testCase.SigFigs, testCase.Records[0].Name, testCase.LongestName)

		err := testOut.Close()
		exception.PanicOnErr(err)
		if fileio.AreEqual(testCase.Out, testCase.Expected) {
			t.Errorf("Error: in in pFasta. printOneSetLines test not as expected.")
		}
	}
}
