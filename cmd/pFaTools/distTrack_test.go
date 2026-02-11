package main

import (
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/wig"
	"os"
	"testing"
)

var distTrackTests = []struct {
	InFile1       string
    InFile2       string
	OutFile      string
    OutName      string
	ChromSizes   string
	ExpectedFile string
	DefaultValue float64
}{
	{InFile1: "testdata/test_sample_input.pfa",
        InFile2: "testdata/test_sample_input_2.pfa"
		OutFile:      "testdata/distTrack.wig",
		ChromSizes:   "testdata/test_sample_input.chrom.sizes",
		ExpectedFile: "testdata/expected.distTrack.wig",
		DefaultValue: -1,
	},
}

func TestDistTrack(t *testing.T) {
	var err error
	var alpha, beta map[string]wig.Wig
	var s DistTrackSettings
	for _, testCase := range distTrackTests {
		s = DistTrackSettings{
			InFile:       testCase.InFile,
            InFile2:       testCase.InFile2,
			OutFile:      testCase.OutFile,
			DefaultValue: testCase.DefaultValue,
            DefaultName: testCase.DefaultName,
		}
		pFaDistTrack(s)
		// TODO: need to create the second input file
		alpha = wig.Read(testCase.OutFile, testCase.ChromSizes, testCase.DefaultValue)
		beta = wig.Read(testCase.ExpectedFile, testCase.ChromSizes, testCase.DefaultValue)
		if !wig.AllEqual(alpha, beta, 1e-6) {
			t.Errorf("Error: DistTrack output was not as expected.\n")
		} else {
			err = os.Remove(testCase.OutFile)
			exception.PanicOnErr(err)
		}
	}
}
