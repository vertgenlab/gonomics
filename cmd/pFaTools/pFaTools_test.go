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
		Chrom:        "chrom1",
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
	Chrom        string
	Region       string
	OutName      string
	OutFile      string
	ExpectedFile string
	Precision    float32
}{
	{InFile: "testdata/test_extract_input.pfa",
		Chrom:        "chrom1",
		Region:       "testdata/test_extractbed_input_region.bed",
		OutName:      "testChr1",
		OutFile:      "testdata/test_extractbed_bed_output.pfa",
		ExpectedFile: "testdata/test_extractbed_expected.pfa",
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
			OutName: testCase.OutName,
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
	OutDir       string
	NumSamples   int
	ExpectedFile string
	SetSeed      uint64
}{
	{InFile: "testdata/test_sample_input.pfa",
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
