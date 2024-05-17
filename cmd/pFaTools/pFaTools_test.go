package main

import (
	"github.com/vertgenlab/gonomics/dna/pDna"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fasta/pFasta"
	"github.com/vertgenlab/gonomics/fileio"
	"os"
	"testing"
	"os"
	"testing"
)

var extractTests = []struct {
	InFile       string
	Chrom        string
	OutName      string
	OutFile      string
	Start        int
	End          int
	ExpectedFile string
	Precision    float32
}{
	{InFile: "testdata/test_extract_input.pfa",
		Chrom:        "chr1",
		OutName:      "testChr1",
		OutFile:      "testdata/test_extract_output.pfa",
		Start:        1,
		End:          3,
		ExpectedFile: "testdata/test_extract_expected.pfa",
		Precision:    1e-3,
	},
}

func TestExtract(t *testing.T) {
	var err error
	var s ExtractSettings
	for _, testCase := range extractTests {
		s = ExtractSettings{
			InFile:  testCase.InFile,
			Chrom:   testCase.Chrom,
			OutName: testCase.OutName,
			OutFile: testCase.OutFile,
			Start:   testCase.Start,
			End:     testCase.End,
		}
		pFaExtract(s)
		observed := pFasta.Read(testCase.OutFile)
		expected := pFasta.Read(testCase.ExpectedFile)
		if !pFasta.AllAreEqual(observed, expected, testCase.Precision) {
			t.Errorf("Error: pFaExtract outFile is not as expected.")
		} else {
			err = os.Remove(testCase.OutFile)
			exception.PanicOnErr(err)
		}
	}
}

var extractBedTests = []struct {
	InFile       string
	Region       string
	TakeCoords   bool
	OutFile      string
	ExpectedFile string
	Precision    float32
}{
	{InFile: "testdata/test_extractbed_input_0.pfa",
		Region:       "testdata/test_extractbed_input_region_0.bed",
		OutFile:      "testdata/test_extractbed_bed_output_0.pfa",
		TakeCoords:   false,
		ExpectedFile: "testdata/test_extractbed_expected_0.pfa",
		Precision:    1e-3,
	},
	{InFile: "testdata/test_extractbed_input_1.pfa",
		Region:       "testdata/test_extractbed_input_region_1.bed",
		OutFile:      "testdata/test_extractbed_bed_output_1.pfa",
		TakeCoords:   true,
		ExpectedFile: "testdata/test_extractbed_expected_1.pfa",
		Precision:    1e-3,
	},
}

func TestExtractBed(t *testing.T) {
	var err error
	var s ExtractBedSettings
	for _, testCase := range extractBedTests {
		s = ExtractBedSettings{
			InFile:  testCase.InFile,
			Region:  testCase.Region,
			TakeCoords: testCase.TakeCoords,
			OutFile: testCase.OutFile,
		}
		pFaExtractBed(s)
		observed := pFasta.Read(testCase.OutFile)
		expected := pFasta.Read(testCase.ExpectedFile)
		if !pFasta.AllAreEqual(observed, expected, testCase.Precision) {
			t.Errorf("Error: pFaExtract outFile is not as expected.")
		} else {
			err = os.Remove(testCase.OutFile)
			exception.PanicOnErr(err)
		}
	}
}

var sampleTests = []struct {
	InFile       string
	Chrom		 string
	OutDir       string
	NumSamples   int
	ExpectedFile string
	SetSeed      int64
}{
	{InFile: "testdata/test_sample_input.pfa",
		Chrom: 		  "chr1",
		OutDir:       "testdata",
		NumSamples:   1,
		ExpectedFile: "testdata/test_sample_expected.fa",
		SetSeed:      7,
	},
}

func TestSample(t *testing.T) {
	var err error
	var s SampleSettings
	for _, testCase := range sampleTests {
		s = SampleSettings{
			InFile:     testCase.InFile,
			Chrom:		testCase.Chrom,
			OutDir:     testCase.OutDir,
			NumSamples: testCase.NumSamples,
			SetSeed:    testCase.SetSeed,
		}
		pFaSample(s)
		for currSample := 0; currSample < s.NumSamples; currSample++ {
			var outName = fmt.Sprintf("%s/sample_%v.fa", s.OutDir, currSample)
			if !fileio.AreEqual(outName, testCase.ExpectedFile) {
				t.Errorf("Error: pFaExtract outFile is not as expected.")
			} else {
				err = os.Remove(outName)
				exception.PanicOnErr(err)
			}
		}
	}
}

var visualizeTests = []struct {
	InFile           string
	OutDir           string
	Start            int
	End              int
	SigFigs          int
	DecimalPlaces    int
	LineLength       int
	Chrom            string
	StartOfAlignment bool
	EndOfAlignment   bool
	ExpectedFile     string
}{
	{InFile: "testdata/test_visualize_input_1.pfa",
		OutDir:           "testdata/test_visualize_output_default.txt",
		Start:            0,
		End:              15,
		SigFigs:          0,
		DecimalPlaces:    5,
		LineLength:       50,
		Chrom:            "chr1butrllllllylong",
		StartOfAlignment: false,
		EndOfAlignment:   false,
		ExpectedFile:     "testdata/test_visualize_expected_default.txt",
	},
	{InFile: "testdata/test_visualize_input_1.pfa",
		OutDir:           "testdata/test_visualize_output_1.txt",
		Start:            0,
		End:              -1,
		SigFigs:          0,
		DecimalPlaces:    15,
		LineLength:       10,
		Chrom:            "chr1butrllllllylong",
		StartOfAlignment: true,
		EndOfAlignment:   true,
		ExpectedFile:     "testdata/test_visualize_expected_1.txt",
	},
	{InFile: "testdata/test_visualize_input_1.pfa",
		OutDir:           "testdata/test_visualize_output_2.txt",
		Start:            0,
		End:              20,
		SigFigs:          0,
		DecimalPlaces:    4,
		LineLength:       7,
		Chrom:            "chr1butrllllllylong",
		StartOfAlignment: true,
		EndOfAlignment:   false,
		ExpectedFile:     "testdata/test_visualize_expected_2.txt",
	},
	{InFile: "testdata/test_visualize_input_2.pfa",
		OutDir:           "testdata/test_visualize_output_3.txt",
		Start:            0,
		End:              -1,
		SigFigs:          10,
		DecimalPlaces:    3,
		LineLength:       50,
		Chrom:            "",
		StartOfAlignment: false,
		EndOfAlignment:   true,
		ExpectedFile:     "testdata/test_visualize_expected_3.txt",
	},
	{InFile: "testdata/test_visualize_normalized_input_1.pfa",
		OutDir:           "testdata/test_visualize_normalized_output_1.txt",
		Start:            2,
		End:              15,
		SigFigs:          0,
		DecimalPlaces:    1,
		LineLength:       6,
		Chrom:            "chr1",
		StartOfAlignment: false,
		EndOfAlignment:   false,
		ExpectedFile:     "testdata/test_visualize_normalized_expected_1.txt",
	},
}

func TestVisualize(t *testing.T) {
	var err error
	var s ExtractSettings
	for _, testCase := range extractTests {
		s = ExtractSettings{
			InFile:  testCase.InFile,
			Chrom:   testCase.Chrom,
			OutName: testCase.OutName,
			OutFile: testCase.OutFile,
			Start:   testCase.Start,
			End:     testCase.End,
		}
		pFaExtract(s)
		observed := pFasta.Read(testCase.OutFile)
		expected := pFasta.Read(testCase.ExpectedFile)
		if !pFasta.AllAreEqual(observed, expected, testCase.Precision) {
			t.Errorf("Error: pFaExtract outFile is not as expected.")
		} else {
			err = os.Remove(testCase.OutFile)
			exception.PanicOnErr(err)
		}
	}
}

var extractBedTests = []struct {
	InFile       string
	Region       string
	TakeCoords   bool
	OutFile      string
	ExpectedFile string
	Precision    float32
}{
	{InFile: "testdata/test_extractbed_input_0.pfa",
		Region:       "testdata/test_extractbed_input_region_0.bed",
		OutFile:      "testdata/test_extractbed_bed_output_0.pfa",
		TakeCoords:   false,
		ExpectedFile: "testdata/test_extractbed_expected_0.pfa",
		Precision:    1e-3,
	},
	{InFile: "testdata/test_extractbed_input_1.pfa",
		Region:       "testdata/test_extractbed_input_region_1.bed",
		OutFile:      "testdata/test_extractbed_bed_output_1.pfa",
		TakeCoords:   true,
		ExpectedFile: "testdata/test_extractbed_expected_1.pfa",
		Precision:    1e-3,
	},
}

func TestExtractBed(t *testing.T) {
	var err error
	var s ExtractBedSettings
	for _, testCase := range extractBedTests {
		s = ExtractBedSettings{
			InFile:     testCase.InFile,
			Region:     testCase.Region,
			TakeCoords: testCase.TakeCoords,
			OutFile:    testCase.OutFile,
		}
		pFaExtractBed(s)
		observed := pFasta.Read(testCase.OutFile)
		expected := pFasta.Read(testCase.ExpectedFile)
		if !pFasta.AllAreEqual(observed, expected, testCase.Precision) {
			t.Errorf("Error: pFaExtract outFile is not as expected.")
		} else {
			err = os.Remove(testCase.OutFile)
			exception.PanicOnErr(err)
		}
	}
}

var sampleTests = []struct {
	InFile       string
	Chrom        string
	OutDir       string
	NumSamples   int
	ExpectedFile string
	SetSeed      int64
}{
	{InFile: "testdata/test_sample_input.pfa",
		Chrom:        "chr1",
		OutDir:       "testdata",
		NumSamples:   1,
		ExpectedFile: "testdata/test_sample_expected.fa",
		SetSeed:      7,
	},
}

func TestSample(t *testing.T) {
	var err error
	var s SampleSettings
	for _, testCase := range sampleTests {
		s = SampleSettings{
			InFile:     testCase.InFile,
			Chrom:      testCase.Chrom,
			OutDir:     testCase.OutDir,
			NumSamples: testCase.NumSamples,
			SetSeed:    testCase.SetSeed,
		}
		pFaSample(s)
		for currSample := 0; currSample < s.NumSamples; currSample++ {
			var outName = fmt.Sprintf("%s/sample_%v.fa", s.OutDir, currSample)
			if !fileio.AreEqual(outName, testCase.ExpectedFile) {
				t.Errorf("Error: pFaExtract outFile is not as expected.")
			} else {
				err = os.Remove(outName)
				exception.PanicOnErr(err)
			}
		}
	}
}

var visualizeTests = []struct {
	InFile           string
	OutDir           string
	Start            int
	End              int
	SigFigs          int
	DecimalPlaces    int
	LineLength       int
	Chrom            string
var visualizeRTests = []struct {
	InFile           string
	OutDir           string
	Start            int
	End              int
	SigFigs          int
	DecimalPlaces    int
	LineLength       int
	Chrom            string
	StartOfAlignment bool
	EndOfAlignment   bool
	ExpectedFile     string
	EndOfAlignment   bool
	ExpectedFile     string
}{
	{InFile: "testdata/test_visualize_input_1.pfa",
		OutDir:           "testdata/test_visualize_output_default.txt",
		Start:            0,
		End:              15,
		SigFigs:          0,
		DecimalPlaces:    5,
		LineLength:       50,
		Chrom:            "chr1butrllllllylong",
		StartOfAlignment: false,
		EndOfAlignment:   false,
		ExpectedFile:     "testdata/test_visualize_expected_default.txt",
		OutDir:           "testdata/test_visualize_output_default.txt",
		Start:            0,
		End:              15,
		SigFigs:          0,
		DecimalPlaces:    5,
		LineLength:       50,
		Chrom:            "chr1butrllllllylong",
		StartOfAlignment: false,
		EndOfAlignment:   false,
		ExpectedFile:     "testdata/test_visualize_expected_default.txt",
	},
	{InFile: "testdata/test_visualize_input_1.pfa",
		OutDir:           "testdata/test_visualize_output_1.txt",
		Start:            0,
		End:              -1,
		SigFigs:          0,
		DecimalPlaces:    15,
		LineLength:       10,
		Chrom:            "chr1butrllllllylong",
		StartOfAlignment: true,
		EndOfAlignment:   true,
		ExpectedFile:     "testdata/test_visualize_expected_1.txt",
		OutDir:           "testdata/test_visualize_output_1.txt",
		Start:            0,
		End:              -1,
		SigFigs:          0,
		DecimalPlaces:    15,
		LineLength:       10,
		Chrom:            "chr1butrllllllylong",
		StartOfAlignment: true,
		EndOfAlignment:   true,
		ExpectedFile:     "testdata/test_visualize_expected_1.txt",
	},
	{InFile: "testdata/test_visualize_input_1.pfa",
		OutDir:           "testdata/test_visualize_output_2.txt",
		Start:            0,
		End:              20,
		SigFigs:          0,
		DecimalPlaces:    4,
		LineLength:       7,
		Chrom:            "chr1butrllllllylong",
		StartOfAlignment: true,
		EndOfAlignment:   false,
		ExpectedFile:     "testdata/test_visualize_expected_2.txt",
		OutDir:           "testdata/test_visualize_output_2.txt",
		Start:            0,
		End:              20,
		SigFigs:          0,
		DecimalPlaces:    4,
		LineLength:       7,
		Chrom:            "chr1butrllllllylong",
		StartOfAlignment: true,
		EndOfAlignment:   false,
		ExpectedFile:     "testdata/test_visualize_expected_2.txt",
	},
	{InFile: "testdata/test_visualize_input_2.pfa",
		OutDir:           "testdata/test_visualize_output_3.txt",
		Start:            0,
		End:              -1,
		SigFigs:          10,
		DecimalPlaces:    3,
		LineLength:       50,
		Chrom:            "",
		StartOfAlignment: false,
		EndOfAlignment:   true,
		ExpectedFile:     "testdata/test_visualize_expected_3.txt",
		OutDir:           "testdata/test_visualize_output_3.txt",
		Start:            0,
		End:              -1,
		SigFigs:          10,
		DecimalPlaces:    3,
		LineLength:       50,
		Chrom:            "",
		StartOfAlignment: false,
		EndOfAlignment:   true,
		ExpectedFile:     "testdata/test_visualize_expected_3.txt",
	},
	{InFile: "testdata/test_visualize_normalized_input_1.pfa",
		OutDir:           "testdata/test_visualize_normalized_output_1.txt",
		Start:            2,
		End:              15,
		SigFigs:          0,
		DecimalPlaces:    1,
		LineLength:       6,
		Chrom:            "chr1",
		StartOfAlignment: false,
		EndOfAlignment:   false,
		ExpectedFile:     "testdata/test_visualize_normalized_expected_1.txt",
		OutDir:           "testdata/test_visualize_normalized_output_1.txt",
		Start:            2,
		End:              15,
		SigFigs:          0,
		DecimalPlaces:    1,
		LineLength:       6,
		Chrom:            "chr1",
		StartOfAlignment: false,
		EndOfAlignment:   false,
		ExpectedFile:     "testdata/test_visualize_normalized_expected_1.txt",
	},
}

func TestVisualize(t *testing.T) {
	var err error
	var s VisualizeSettings

	for _, testCase := range visualizeTests {
		s = VisualizeSettings{
			InFile:           testCase.InFile,
			OutDir:           testCase.OutDir,
			Start:            testCase.Start,
			End:              testCase.End,
			SigFigs:          testCase.SigFigs,
			DecimalPlaces:    testCase.DecimalPlaces,
			LineLength:       testCase.LineLength,
			Chrom:            testCase.Chrom,
			StartOfAlignment: testCase.StartOfAlignment,
			EndOfAlignment:   testCase.EndOfAlignment,
		}
		pFaVisualize(s)


		if !fileio.AreEqual(testCase.OutDir, testCase.ExpectedFile) {
			t.Errorf("Error: pFaVisualise output not as expected.")
		} else {
			err = os.Remove(testCase.OutDir)
			exception.PanicOnErr(err)
		}
	}
}
