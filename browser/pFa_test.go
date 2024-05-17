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
	Records: []pFasta.PFasta{{Name: "chr1",
		Seq: []pDna.Float32Base{
			{
				A: 0.23857,
				C: 0.3323,
				G: 0.44958,
				T: 0.139448,
			},
			{
				A: 0.334,
				C: 0.44239,
				G: 0.134875,
				T: 0.12394,
			},
			{
				A: 0.384398,
				C: 0.59723,
				G: 0.3289,
				T: 0.8325,
			},
			{
				A: 0.488532,
				C: 0.23858,
				G: 0.548523,
				T: 0.5473,
			},
			{
				A: 0.92323,
				C: 0.237,
				G: 0.1747,
				T: 0.5839,
			},
			{
				A: 0.483284,
				C: 0.123,
				G: 0.25388,
				T: 0.98243,
			},
			{
				A: 0.00834,
				C: 0.5288,
				G: 0.58001,
				T: 0.4892,
			},
			{
				A: 0.5688,
				C: 0.278,
				G: 0.1937,
				T: 0.1111,
			},
			{
				A: 0.42397,
				C: 0.2358,
				G: 0.984,
				T: 0.3823,
			},
			{
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
		} else {
			fileio.EasyRemove(testCase.Out)
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
		} else {
			fileio.EasyRemove(testCase.OutFile)
		}
	}
}

var PFaVisualizerTests = []struct {
	InFile           string
	OutFile          string
	Start            int
	End              int
	StartOfAlignment bool
	EndOfAlignment   bool
	SigFigs          int
	DecimalPlaces    int
	LineLength       int
	SeqName          string
	Expected         string
}{{InFile: "testdata/pfa_PFaVisualiser_input_toy_1.pfa",
	OutFile:          "testdata/pfa_PFaVisualiser_output_toy_1.txt",
	Start:            4,
	End:              21,
	StartOfAlignment: false,
	EndOfAlignment:   false,
	SigFigs:          4,
	DecimalPlaces:    7,
	LineLength:       5,
	SeqName:          "chr1",
	Expected:         "testdata/pfa_PFaVisualiser_expected_toy_1.txt",
},
	{InFile: "testdata/pfa_PFaVisualiser_input_toy_1.pfa",
		OutFile:          "testdata/pfa_PFaVisualiser_output_toy_2.txt",
		Start:            4,
		End:              21,
		StartOfAlignment: false,
		EndOfAlignment:   false,
		SigFigs:          0,
		DecimalPlaces:    7,
		LineLength:       5,
		SeqName:          "chr1",
		Expected:         "testdata/pfa_PFaVisualiser_expected_toy_2.txt",
	},
	{InFile: "testdata/pfa_PFaVisualiser_input_toy_1.pfa",
		OutFile:          "testdata/pfa_PFaVisualiser_output_toy_3.txt",
		Start:            4,
		End:              21,
		StartOfAlignment: false,
		EndOfAlignment:   false,
		SigFigs:          0,
		DecimalPlaces:    4,
		LineLength:       5,
		SeqName:          "chr1",
		Expected:         "testdata/pfa_PFaVisualiser_expected_toy_3.txt",
	},
	{InFile: "testdata/pfa_PFaVisualiser_normalised_input_toy_1.pfa",
		OutFile:          "testdata/pfa_PFaVisualiser_normalised_output_toy_1.txt",
		Start:            6,
		End:              13,
		StartOfAlignment: false,
		EndOfAlignment:   false,
		SigFigs:          2,
		DecimalPlaces:    5,
		LineLength:       4,
		SeqName:          "chr1",
		Expected:         "testdata/pfa_PFaVisualiser_normalised_expected_toy_1.txt",
	},
	{InFile: "testdata/pfa_PFaVisualiser_normalised_input_toy_1.pfa",
		OutFile:          "testdata/pfa_PFaVisualiser_normalised_output_toy_2.txt",
		Start:            6,
		End:              13,
		StartOfAlignment: false,
		EndOfAlignment:   false,
		SigFigs:          0,
		DecimalPlaces:    5,
		LineLength:       4,
		SeqName:          "chr1",
		Expected:         "testdata/pfa_PFaVisualiser_normalised_expected_toy_2.txt",
	},
	{InFile: "testdata/pfa_PFaVisualiser_normalised_input_toy_2.pfa",
		OutFile:          "testdata/pfa_PFaVisualiser_normalised_output_toy_3.txt",
		Start:            6,
		End:              13,
		StartOfAlignment: false,
		EndOfAlignment:   false,
		SigFigs:          2,
		DecimalPlaces:    5,
		LineLength:       4,
		SeqName:          "chr1butlikerealllllllllllllylong",
		Expected:         "testdata/pfa_PFaVisualiser_normalised_expected_toy_3.txt",
	},
}

func TestPFaVisualizer(t *testing.T) {
	for _, testCase := range PFaVisualizerTests {
		PFaVisualizer(testCase.InFile, testCase.OutFile, testCase.Start, testCase.End, testCase.StartOfAlignment, testCase.EndOfAlignment, testCase.SigFigs, testCase.DecimalPlaces, testCase.LineLength, testCase.SeqName)
		if !fileio.AreEqual(testCase.OutFile, testCase.Expected) {
			t.Errorf("Error: in browser. PFaVisualiser test not as expected.")
		} else {
			fileio.EasyRemove(testCase.OutFile)
		}

	}
}

var PFaVisualizerRTests = []struct {
	InFile           string
	OutFile          string
	Start            int
	End              int
	StartOfAlignment bool
	EndOfAlignment   bool
	SigFigs          int
	DecimalPlaces    int
	LineLength       int
	SeqName          string
	Expected         string
}{
	{InFile: "testdata/pfa_PFaVisualiser_input_toy_1.pfa",
		OutFile:          "testdata/pfa_PFaVisualiserR_output_toy_1.txt",
		Start:            4,
		End:              21,
		StartOfAlignment: false,
		EndOfAlignment:   false,
		SigFigs:          4,
		DecimalPlaces:    7,
		LineLength:       5,
		SeqName:          "chr1",
		Expected:         "testdata/pfa_PFaVisualiserR_expected_toy_1.txt",
	},
	{InFile: "testdata/pfa_PFaVisualiser_input_toy_1.pfa",
		OutFile:          "testdata/pfa_PFaVisualiserR_output_toy_2.txt",
		Start:            4,
		End:              21,
		StartOfAlignment: false,
		EndOfAlignment:   false,
		SigFigs:          0,
		DecimalPlaces:    7,
		LineLength:       5,
		SeqName:          "chr1",
		Expected:         "testdata/pfa_PFaVisualiserR_expected_toy_2.txt",
	},
	{InFile: "testdata/pfa_PFaVisualiser_input_toy_1.pfa",
		OutFile:          "testdata/pfa_PFaVisualiserR_output_toy_3.txt",
		Start:            4,
		End:              21,
		StartOfAlignment: false,
		EndOfAlignment:   false,
		SigFigs:          0,
		DecimalPlaces:    4,
		LineLength:       5,
		SeqName:          "chr1",
		Expected:         "testdata/pfa_PFaVisualiserR_expected_toy_3.txt",
	},
	{InFile: "testdata/pfa_PFaVisualiser_normalised_input_toy_1.pfa",
		OutFile:          "testdata/pfa_PFaVisualiserR_normalised_output_toy_1.txt",
		Start:            6,
		End:              13,
		StartOfAlignment: false,
		EndOfAlignment:   false,
		SigFigs:          2,
		DecimalPlaces:    5,
		LineLength:       4,
		SeqName:          "chr1",
		Expected:         "testdata/pfa_PFaVisualiserR_normalised_expected_toy_1.txt",
	},
	{InFile: "testdata/pfa_PFaVisualiser_normalised_input_toy_1.pfa",
		OutFile:          "testdata/pfa_PFaVisualiserR_normalised_output_toy_2.txt",
		Start:            6,
		End:              13,
		StartOfAlignment: false,
		EndOfAlignment:   false,
		SigFigs:          0,
		DecimalPlaces:    5,
		LineLength:       4,
		SeqName:          "chr1",
		Expected:         "testdata/pfa_PFaVisualiserR_normalised_expected_toy_2.txt",
	},
	{InFile: "testdata/pfa_PFaVisualiser_normalised_input_toy_2.pfa",
		OutFile:          "testdata/pfa_PFaVisualiserR_normalised_output_toy_3.txt",
		Start:            6,
		End:              13,
		StartOfAlignment: false,
		EndOfAlignment:   false,
		SigFigs:          2,
		DecimalPlaces:    5,
		LineLength:       4,
		SeqName:          "chr1butlikerealllllllllllllylong",
		Expected:         "testdata/pfa_PFaVisualiserR_normalised_expected_toy_3.txt",
	},
}

func TestPFaVisualizerR(t *testing.T) {
	for _, testCase := range PFaVisualizerRTests {
		PFaVisualizerR(testCase.InFile, testCase.OutFile, testCase.Start, testCase.End, testCase.StartOfAlignment, testCase.EndOfAlignment, testCase.SigFigs, testCase.DecimalPlaces, testCase.LineLength, testCase.SeqName)
		if !fileio.AreEqual(testCase.OutFile, testCase.Expected) {
			t.Errorf("Error: in browser. PFaVisualiserR test not as expected.")
		} else {
			fileio.EasyRemove(testCase.OutFile)
		}
	}
}
