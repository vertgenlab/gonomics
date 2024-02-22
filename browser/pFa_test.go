
package browser

import (
	"testing"
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
	LineLength     int
	SetOfLinesIdx  int
	NumIters int
	LineA	[]float32
	LineC	[]float32
	LineG	[]float32
	LineT	[]float32
	Start int
	Records []PFasta
	Out	*fileio.EasyWriter
	SigFigs int
	Expected	string
}{{LineLength: 30,
	SetOfLinesIdx: 0,
	NumIters: 30,
	LineA: make([]float32, 30),
	LineC: make([]float32, 30),
	LineG: make([]float32, 30),
	LineT: make([]float32, 30),
	Start: 2,
	Records: ,
	Out: fileio.EasyCreate('testdata_tools/printONeSetLinesTest_1.txt'),
	SigFigs: 3,
	Expected: 'testdata_tools/printONeSetLinesTest_1_Expected.txt',
	}
}

// printOneSetLines(lineLength int, setOfLinesIdx int, numIters int, lineA []float32, lineC []float32, lineG []float32, lineT []float32, start int, records []PFasta, out *fileio.EasyWriter, sigFigs int) {

func TestPrintOneSetLines(t *testing.T) {
	for _, testCase := range PrintOneSetLinesTests {
		printOneSetLines(testCase.lineLength, testCase.setOfLinesIdx, testCase.numIters, testCase.lineA, testCase.lineC, testCase.lineG, testCase.lineT, testCase.start, testCase.records, testCase.out, testCase.sigFigs)

		if resA != testCase.Expected[0] || resC != testCase.Expected[1] || resG != testCase.Expected[2] || resT != testCase.Expected[3] {
			fmt.Printf("%v, %v, %v, %v\n", resA, resC, resG, resT)

			t.Errorf("Error: in in pFasta. printOneSetLines test not as expected.")
		}
	}
}
