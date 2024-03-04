package browser

import (
	"github.com/vertgenlab/gonomics/dna/pDna"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fasta/pFasta"
	"github.com/vertgenlab/gonomics/fileio"
	"testing"
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
			t.Errorf("Error: in browser. getBaseProbsAtPos test not as expected.")
		}
	}

	pFaInput := "testdata/pfa_PFaVisualiser_input_toy_1.pfa"
	Records := []pFasta.PFasta{
		pFasta.PFasta{Name: "chr2",
			Seq: []pDna.Float32Base{
				pDna.Float32Base{
					A: 0.1,
					C: 0.2,
					G: 0.3,
					T: 0.4}}},
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
					A: 0.32984,
					C: 0.123,
					G: 0.482,
					T: 0.1,
				},
				pDna.Float32Base{
					A: 0.2,
					C: 0.3495,
					G: 0.43012,
					T: 0.239,
				},
				pDna.Float32Base{
					A: 0.2349,
					C: 0.91273,
					G: 0.91237,
					T: 0.127,
				},
				pDna.Float32Base{
					A: 0.12398,
					C: 0.9438,
					G: 0.1238,
					T: 0.93458,
				},
				pDna.Float32Base{
					A: 0.384,
					C: 0.1874,
					G: 0.23498,
					T: 0.32489,
				},
				pDna.Float32Base{
					A: 0.23982,
					C: 0.4571,
					G: 0.349,
					T: 0.982,
				},
				pDna.Float32Base{
					A: 0.23948,
					C: 0.4557,
					G: 0.65832,
					T: 0.458,
				},
				pDna.Float32Base{
					A: 0.85,
					C: 0.5703,
					G: 0.2047,
					T: 0.29047,
				},
				pDna.Float32Base{
					A: 0.0384,
					C: 0.03874,
					G: 0.9273,
					T: 0.0837,
				},
				pDna.Float32Base{
					A: 0.0582,
					C: 0.31495,
					G: 0.617,
					T: 0.791,
				},
				pDna.Float32Base{
					A: 0.9813,
					C: 0.1239,
					G: 0.348,
					T: 0.29438,
				},
				pDna.Float32Base{
					A: 0.34280,
					C: 0.9872,
					G: 0.298,
					T: 0.4789,
				},
				pDna.Float32Base{
					A: 0.4892,
					C: 0.48925,
					G: 0.190,
					T: 0.549,
				},
				pDna.Float32Base{
					A: 0.3480,
					C: 0.00124,
					G: 0.60693,
					T: 0.09912,
				},
				pDna.Float32Base{
					A: 0.8327,
					C: 0.0459,
					G: 0.384,
					T: 0.4983,
				},
				pDna.Float32Base{
					A: 0.5938,
					C: 0.9481,
					G: 0.1925,
					T: 0.9871,
				},
				pDna.Float32Base{
					A: 0.00123,
					C: 0.0104,
					G: 0.6902,
					T: 0.2398,
				},
				pDna.Float32Base{
					A: 0.0383,
					C: 0.0948,
					G: 0.028,
					T: 0.981,
				}, pDna.Float32Base{
					A: 0.763,
					C: 0.304,
					G: 0.411,
					T: 0.1211,
				},
				pDna.Float32Base{
					A: 0.884,
					C: 0.3955,
					G: 0.8387,
					T: 0.3884,
				},
				pDna.Float32Base{
					A: 0.042387,
					C: 0.587,
					G: 0.32871,
					T: 0.73224,
				},
			},
		},
	}
	pFasta.Write(pFaInput, Records)
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
	DecimalPlaces int
	LongestName   int
	Expected      string
}{{LineLength: 10,
	SetOfLinesIdx: 0,
	NumIters:      6,
	LineA:         make([]float32, 30),
	LineC:         make([]float32, 30),
	LineG:         make([]float32, 30),
	LineT:         make([]float32, 30),
	Start:         2,
	Records: []pFasta.PFasta{pFasta.PFasta{Name: "chr1",
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
	Out:           "testdata/pfa_printOneSetLines_output_toy_1.txt",
	SigFigs:       3,
	DecimalPlaces: 7,
	LongestName:   15,
	Expected:      "testdata/pfa_printOneSetLines_expected_toy_1.txt",
},
}

func TestPrintOneSetLines(t *testing.T) {
	for _, testCase := range PrintOneSetLinesTests {
		testOut := fileio.EasyCreate(testCase.Out)

		printOneSetLines(testCase.LineLength, testCase.SetOfLinesIdx, testCase.NumIters,
			testCase.LineA, testCase.LineC, testCase.LineG, testCase.LineT, testCase.Start,
			testCase.Records[0], testOut, testCase.SigFigs, testCase.DecimalPlaces, testCase.LongestName)

		err := testOut.Close()
		exception.PanicOnErr(err)
		if !fileio.AreEqual(testCase.Out, testCase.Expected) {
			t.Errorf("Error: in browser. printOneSetLines test not as expected.")
		}
	}
}

var PrintAllSetsTests = []struct {
	InFile        string
	OutFile       string
	Start         int
	End           int
	LineLength    int
	SigFigs       int
	DecimalPlaces int
	PFaIdx        int
	Expected      string
}{
	{InFile: "testdata/pfa_PFaVisualiser_input_toy_1.pfa",
		OutFile:       "testdata/pfa_printAllSets_output_toy_1.txt",
		Start:         4,
		End:           21,
		LineLength:    5,
		SigFigs:       4,
		DecimalPlaces: 7,
		PFaIdx:        1,
		Expected:      "testdata/pfa_printAllSets_expected_toy_1.txt",
	},
}

func TestPrintAllSets(t *testing.T) {
	for _, testCase := range PrintAllSetsTests {
		records := pFasta.Read(testCase.InFile)
		out := fileio.EasyCreate(testCase.OutFile)
		var err error
		printAllSets(out, err, records[testCase.PFaIdx], testCase.Start, testCase.End, testCase.LineLength, testCase.SigFigs, testCase.DecimalPlaces)
		err = out.Close()
		exception.PanicOnErr(err)
		if !fileio.AreEqual(testCase.OutFile, testCase.Expected) {
			t.Errorf("Error: in browser. PrintAllSets test not as expected.")
		}
	}
}

var PFaVisualizerTests = []struct {
	InFile         string
	OutFile        string
	Start          int
	End            int
	EndOfAlignment bool
	SigFigs        int
	DecimalPlaces  int
	LineLength     int
	SeqName        string
	Expected       string
}{{InFile: "testdata/pfa_PFavisualiser_input_toy_1.pfa",
	OutFile:        "testdata/pfa_PFaVisualiser_output_toy_1.txt",
	EndOfAlignment: false,
	Start:          4,
	End:            21,
	SigFigs:        4,
	DecimalPlaces:  7,
	LineLength:     5,
	SeqName:        "chr1",
	Expected:       "testdata/pfa_PFaVisualiser_expected_toy_1.txt",
},
	{InFile: "testdata/pfa_PFaVisualiser_input_toy_1.pfa",
		OutFile:        "testdata/pfa_PFaVisualiser_output_toy_2.txt",
		EndOfAlignment: false,
		Start:          4,
		End:            21,
		SigFigs:        0,
		DecimalPlaces:  7,
		LineLength:     5,
		SeqName:        "chr1",
		Expected:       "testdata/pfa_PFaVisualiser_expected_toy_2.txt",
	},
	{InFile: "testdata/pfa_PFaVisualiser_input_toy_1.pfa",
		OutFile:        "testdata/pfa_PFaVisualiser_output_toy_3.txt",
		EndOfAlignment: false,
		Start:          4,
		End:            21,
		SigFigs:        0,
		DecimalPlaces:  4,
		LineLength:     5,
		SeqName:        "chr1",
		Expected:       "testdata/pfa_PFaVisualiser_expected_toy_3.txt",
	},
}

func TestPFaVisualizer(t *testing.T) {
	for _, testCase := range PFaVisualizerTests {
		PFaVisualizer(testCase.InFile, testCase.OutFile, testCase.Start, testCase.End, testCase.EndOfAlignment, testCase.SigFigs, testCase.DecimalPlaces, testCase.LineLength, testCase.SeqName)
		if !fileio.AreEqual(testCase.OutFile, testCase.Expected) {
			t.Errorf("Error: in browser. PFaVisualiser test not as expected.")
		}
	}
}
