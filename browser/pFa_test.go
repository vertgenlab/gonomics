package browser

import (
	"github.com/vertgenlab/gonomics/dna/pDna"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fasta/pFasta"
	"github.com/vertgenlab/gonomics/fileio"
	"testing"
)

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

	pFaInputnorm := "testdata/pfa_PFaVisualiser_normalised_input_toy_2.pfa"
	Recordsnorm := []pFasta.PFasta{
		pFasta.PFasta{Name: "chr2",
			Seq: []pDna.Float32Base{
				pDna.Float32Base{
					A: 0.1,
					C: 0.2,
					G: 0.3,
					T: 0.4}}},
		pFasta.PFasta{Name: "chr1butlikerealllllllllllllylong",
			Seq: []pDna.Float32Base{
				pDna.Float32Base{
					A: 0.20564,
					C: 0.2865,
					G: 0.38760,
					T: 0.120224,
				},
				pDna.Float32Base{
					A: 0.323,
					C: 0.42701,
					G: 0.13029,
					T: 0.1197,
				},
				pDna.Float32Base{
					A: 0.179326,
					C: 0.278674,
					G: 0.1535,
					T: 0.3885,
				},
				pDna.Float32Base{
					A: 0.267992,
					C: 0.13088,
					G: 0.300928,
					T: 0.3002,
				},
				pDna.Float32Base{
					A: 0.481059,
					C: 0.124,
					G: 0.091041,
					T: 0.3039,
				},
				pDna.Float32Base{
					A: 0.26229,
					C: 0.0668,
					G: 0.13778,
					T: 0.53313,
				},
				pDna.Float32Base{
					A: 0.00819,
					C: 0.3783,
					G: 0.36011,
					T: 0.2534,
				},
				pDna.Float32Base{
					A: 0.4939,
					C: 0.24142,
					G: 0.1682,
					T: 0.09648,
				},
				pDna.Float32Base{
					A: 0.20926,
					C: 0.1164,
					G: 0.486,
					T: 0.18834,
				},
				pDna.Float32Base{
					A: 0.32963,
					C: 0.11237,
					G: 0.468,
					T: 0.09,
				},
				pDna.Float32Base{
					A: 0.2,
					C: 0.25104,
					G: 0.35296,
					T: 0.196,
				},
				pDna.Float32Base{
					A: 0.1074,
					C: 0.41734,
					G: 0.41726,
					T: 0.058,
				},
				pDna.Float32Base{
					A: 0.05832,
					C: 0.4439,
					G: 0.05822,
					T: 0.43956,
				},
				pDna.Float32Base{
					A: 0.339,
					C: 0.1657,
					G: 0.20771,
					T: 0.28759,
				},
				pDna.Float32Base{
					A: 0.1182634,
					C: 0.225737,
					G: 0.172,
					T: 0.484,
				},
				pDna.Float32Base{
					A: 0.13294,
					C: 0.25065,
					G: 0.36341,
					T: 0.253,
				},
				pDna.Float32Base{
					A: 0.44,
					C: 0.29766,
					G: 0.1069,
					T: 0.15544,
				},
			},
		},
	}
	pFasta.Write(pFaInputnorm, Recordsnorm)
}

var PFaVisualizerTests = []struct {
	InFile         string
	OutFile        string
	Start          int
	End            int
	StartOfAlignment bool
	EndOfAlignment bool
	SigFigs        int
	DecimalPlaces  int
	LineLength     int
	SeqName        string
	Expected       string
}{{InFile: "testdata/pfa_PFaVisualiser_input_toy_1.pfa",
	OutFile:        "testdata/pfa_PFaVisualiser_output_toy_1.txt",
	Start:          4,
	End:            21,
	StartOfAlignment: false,
	EndOfAlignment: false,
	SigFigs:        4,
	DecimalPlaces:  7,
	LineLength:     5,
	SeqName:        "chr1",
	Expected:       "testdata/pfa_PFaVisualiser_expected_toy_1.txt",
},
	{InFile: "testdata/pfa_PFaVisualiser_input_toy_1.pfa",
		OutFile:        "testdata/pfa_PFaVisualiser_output_toy_2.txt",
		Start:          4,
		End:            21,
		StartOfAlignment: false,
		EndOfAlignment: false,
		SigFigs:        0,
		DecimalPlaces:  7,
		LineLength:     5,
		SeqName:        "chr1",
		Expected:       "testdata/pfa_PFaVisualiser_expected_toy_2.txt",
	},
	{InFile: "testdata/pfa_PFaVisualiser_input_toy_1.pfa",
		OutFile:        "testdata/pfa_PFaVisualiser_output_toy_3.txt",
		Start:          4,
		End:            21,
		StartOfAlignment: false,
		EndOfAlignment: false,
		SigFigs:        0,
		DecimalPlaces:  4,
		LineLength:     5,
		SeqName:        "chr1",
		Expected:       "testdata/pfa_PFaVisualiser_expected_toy_3.txt",
	},
	{InFile: "testdata/pfa_PFaVisualiser_normalised_input_toy_1.pfa",
		OutFile:        "testdata/pfa_PFaVisualiser_normalised_output_toy_1.txt",
		Start:          6,
		End:            13,
		StartOfAlignment: false,
		EndOfAlignment: false,
		SigFigs:        2,
		DecimalPlaces:  5,
		LineLength:     4,
		SeqName:        "chr1",
		Expected:       "testdata/pfa_PFaVisualiser_normalised_expected_toy_1.txt",
	},
	{InFile: "testdata/pfa_PFaVisualiser_normalised_input_toy_1.pfa",
		OutFile:        "testdata/pfa_PFaVisualiser_normalised_output_toy_2.txt",
		Start:          6,
		End:            13,
		StartOfAlignment: false,
		EndOfAlignment: false,
		SigFigs:        0,
		DecimalPlaces:  5,
		LineLength:     4,
		SeqName:        "chr1",
		Expected:       "testdata/pfa_PFaVisualiser_normalised_expected_toy_2.txt",
	},
	{InFile: "testdata/pfa_PFaVisualiser_normalised_input_toy_2.pfa",
		OutFile:        "testdata/pfa_PFaVisualiser_normalised_output_toy_3.txt",
		Start:          6,
		End:            13,
		StartOfAlignment: false,
		EndOfAlignment: false,
		SigFigs:        2,
		DecimalPlaces:  5,
		LineLength:     4,
		SeqName:        "chr1butlikerealllllllllllllylong",
		Expected:       "testdata/pfa_PFaVisualiser_normalised_expected_toy_3.txt",
	},
}

func TestPFaVisualizer(t *testing.T) {
	for _, testCase := range PFaVisualizerTests {
		PFaVisualizer(testCase.InFile, testCase.OutFile, testCase.Start, testCase.End, testCase.StartOfAlignment, testCase.EndOfAlignment, testCase.SigFigs, testCase.DecimalPlaces, testCase.LineLength, testCase.SeqName)
		if !fileio.AreEqual(testCase.OutFile, testCase.Expected) {
			t.Errorf("Error: in browser. PFaVisualiser test not as expected.")
		}
	}
}
