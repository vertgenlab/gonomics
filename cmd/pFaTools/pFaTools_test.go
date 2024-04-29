package main

import (
	"fmt"
	"github.com/vertgenlab/gonomics/dna/pDna"
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
	var s VisualizeSettings

	pFaInput1 := "testdata/test_visualize_input_1.pfa"
	Records1 := []pFasta.PFasta{
		{Name: "chr2",
			Seq: []pDna.Float32Base{
				{
					A: 0.1,
					C: 0.2,
					G: 0.3,
					T: 0.4}}},
		{Name: "chr1butrllllllylong",
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
					A: 0.32984,
					C: 0.123,
					G: 0.482,
					T: 0.1,
				},
				{
					A: 0.2,
					C: 0.3495,
					G: 0.43012,
					T: 0.239,
				},
				{
					A: 0.2349,
					C: 0.91273,
					G: 0.91237,
					T: 0.127,
				},
				{
					A: 0.12398,
					C: 0.9438,
					G: 0.1238,
					T: 0.93458,
				},
				{
					A: 0.384,
					C: 0.1874,
					G: 0.23498,
					T: 0.32489,
				},
				{
					A: 0.23982,
					C: 0.4571,
					G: 0.349,
					T: 0.982,
				},
				{
					A: 0.23948,
					C: 0.4557,
					G: 0.65832,
					T: 0.458,
				},
				{
					A: 0.85,
					C: 0.5703,
					G: 0.2047,
					T: 0.29047,
				},
				{
					A: 0.0384,
					C: 0.03874,
					G: 0.9273,
					T: 0.0837,
				},
				{
					A: 0.0582,
					C: 0.31495,
					G: 0.617,
					T: 0.791,
				},
				{
					A: 0.9813,
					C: 0.1239,
					G: 0.348,
					T: 0.29438,
				},
				{
					A: 0.34280,
					C: 0.9872,
					G: 0.298,
					T: 0.4789,
				},
				{
					A: 0.4892,
					C: 0.48925,
					G: 0.190,
					T: 0.549,
				},
				{
					A: 0.3480,
					C: 0.00124,
					G: 0.60693,
					T: 0.09912,
				},
				{
					A: 0.8327,
					C: 0.0459,
					G: 0.384,
					T: 0.4983,
				},
				{
					A: 0.5938,
					C: 0.9481,
					G: 0.1925,
					T: 0.9871,
				},
				{
					A: 0.00123,
					C: 0.0104,
					G: 0.6902,
					T: 0.2398,
				},
				{
					A: 0.0383,
					C: 0.0948,
					G: 0.028,
					T: 0.981,
				}, {
					A: 0.763,
					C: 0.304,
					G: 0.411,
					T: 0.1211,
				},
				{
					A: 0.884,
					C: 0.3955,
					G: 0.8387,
					T: 0.3884,
				},
				{
					A: 0.042387,
					C: 0.587,
					G: 0.32871,
					T: 0.73224,
				},
			},
		},
	}
	pFasta.Write(pFaInput1, Records1)

	pFaInput2 := "testdata/test_visualize_input_2.pfa"
	Records2 := []pFasta.PFasta{
		{Name: "chr1",
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
					A: 0.32984,
					C: 0.123,
					G: 0.482,
					T: 0.1,
				},
				{
					A: 0.2,
					C: 0.3495,
					G: 0.43012,
					T: 0.239,
				},
				{
					A: 0.2349,
					C: 0.91273,
					G: 0.91237,
					T: 0.127,
				},
				{
					A: 0.12398,
					C: 0.9438,
					G: 0.1238,
					T: 0.93458,
				},
				{
					A: 0.384,
					C: 0.1874,
					G: 0.23498,
					T: 0.32489,
				},
				{
					A: 0.23982,
					C: 0.4571,
					G: 0.349,
					T: 0.982,
				},
				{
					A: 0.23948,
					C: 0.4557,
					G: 0.65832,
					T: 0.458,
				},
				{
					A: 0.85,
					C: 0.5703,
					G: 0.2047,
					T: 0.29047,
				},
				{
					A: 0.0384,
					C: 0.03874,
					G: 0.9273,
					T: 0.0837,
				},
				{
					A: 0.0582,
					C: 0.31495,
					G: 0.617,
					T: 0.791,
				},
				{
					A: 0.9813,
					C: 0.1239,
					G: 0.348,
					T: 0.29438,
				},
				{
					A: 0.34280,
					C: 0.9872,
					G: 0.298,
					T: 0.4789,
				},
				{
					A: 0.4892,
					C: 0.48925,
					G: 0.190,
					T: 0.549,
				},
				{
					A: 0.3480,
					C: 0.00124,
					G: 0.60693,
					T: 0.09912,
				},
				{
					A: 0.8327,
					C: 0.0459,
					G: 0.384,
					T: 0.4983,
				},
				{
					A: 0.5938,
					C: 0.9481,
					G: 0.1925,
					T: 0.9871,
				},
				{
					A: 0.00123,
					C: 0.0104,
					G: 0.6902,
					T: 0.2398,
				},
				{
					A: 0.0383,
					C: 0.0948,
					G: 0.028,
					T: 0.981,
				}, {
					A: 0.763,
					C: 0.304,
					G: 0.411,
					T: 0.1211,
				},
				{
					A: 0.884,
					C: 0.3955,
					G: 0.8387,
					T: 0.3884,
				},
				{
					A: 0.042387,
					C: 0.587,
					G: 0.32871,
					T: 0.73224,
				},
			},
		},
	}
	pFasta.Write(pFaInput2, Records2)

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
