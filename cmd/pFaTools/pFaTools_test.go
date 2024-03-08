package main

import (
	"fmt"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fasta/pFasta"
	"github.com/vertgenlab/gonomics/fileio"
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
	InFile	string
	OutDir	string
	Start	int
	End		int
	SigFigs int
	DecimalPlaces int
	LineLength int
	Chrom string
	StartOfAlignment bool
	EndOfAlignment bool
}{
	{InFile: "pfa_PFaVisualiser_input_toy_1.pfa",
	OutDir: "testdata/test_visualise_output_1.txt",
	Start: 4,
	End: 21,
	SigFigs: 0,
	DecimalPlaces: 4,
	LineLength: 5,
	Chrom: "chr1",
	StartOfAlignment: false,
	EndOfAlignment: false,
	Expected: "testdata/pfa_PFaVisualiser_expected_toy_3.txt"
	},
	{InFile: "pfa_PFaVisualiser_input_toy_2.pfa",
	OutDir: "testdata/test_visualise_output_1.txt",
	Start: 4,
	End: -1,
	SigFigs: 2,
	DecimalPlaces: 5,
	LineLength: 4,
	Chrom: "chr1",
	StartOfAlignment: false,
	EndOfAlignment: true,
	Expected: "testdata/pfa_PFaVisualiser_expected_toy_1.txt"
	},
	{InFile: "pfa_PFaVisualiser_normalised_input_toy_1.pfa",
	OutDir: "testdata/test_visualise_output_1.txt",
	Start: 0,
	End: 13,
	SigFigs: 2,
	DecimalPlaces: 5,
	LineLength: 4,
	Chrom: "chr1",
	StartOfAlignment: true,
	EndOfAlignment: false,
	Expected: "testdata/pfa_PFaVisualiser_normalised_expected_toy_1.txt"
	},
	{InFile: "pfa_PFaVisualiser_normalised_input_toy_2.pfa",
	OutDir: "testdata/test_visualise_output_1.txt",
	Start: 0,
	End: -1,
	SigFigs: 2,
	DecimalPlaces: 5,
	LineLength: 4,
	Chrom: "chr1",
	StartOfAlignment: true,
	EndOfAlignment: true,
	Expected: "testdata/pfa_PFaVisualiser_normalised_expected_toy_1.txt"
	},
}

func TestVisualize(t *testing.T) {
	var err error
	var s VisualizeSettings
	for _, testCase := range visualizeTests {
		s = SampleSettings{
			InFile:     testCase.InFile,
			OutDir:     testCase.OutDir,
			Start: testCase.Start,
			End: testCase.End,
			DecimalPlaces: testCase.DecimalPlaces,
			LineLength: testCase.LineLength,
			Chrom:		testCase.Chrom,
			StartOfAlignment: testCase.StartOfAlignment
			EndOfAlignment: testCase.EndOfAlignment
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
