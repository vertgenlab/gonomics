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
	DefaultValue float64
    OutName      string
	ChromSizes   string
	ExpectedFile string
	
}{
	{InFile1: "testdata/test_distTrack_input_A_1.pfa",
        InFile2: "testdata/test_distTrack_input_B_1.pfa",
		OutName: "chr1.dist",
		OutFile:      "testdata/test_distTrack_1_out.wig",
		ChromSizes:   "testdata/test_distTrack_1.chrom.sizes",
		ExpectedFile: "testdata/test_distTrack_1_expected.wig",
		DefaultValue: -1,
	},
	{InFile1: "testdata/test_distTrack_input_A_2.pfa",
        InFile2: "testdata/test_distTrack_input_B_2.pfa",
		OutName: "chr1.dist",
		OutFile:      "testdata/test_distTrack_2_out.wig",
		ChromSizes:   "testdata/test_distTrack_1.chrom.sizes",
		ExpectedFile: "testdata/test_distTrack_2_expected.wig",
		DefaultValue: -1,
	},
	{InFile1: "testdata/test_distTrackFasta_input_A_1.pfa",
        InFile2: "testdata/test_distTrackFasta_input_B_1.fa",
		OutName: "chr1.dist",
		OutFile:      "testdata/test_distTrackFasta_1_out.wig",
		ChromSizes:   "testdata/test_distTrack_1.chrom.sizes",
		ExpectedFile: "testdata/test_distTrackFasta_1_expected.wig",
		DefaultValue: -1,
	},
	{InFile1: "testdata/test_distTrackFasta_input_A_2.pfa",
        InFile2: "testdata/test_distTrackFasta_input_B_2.fa",
		OutName: "chr1.dist",
		OutFile:      "testdata/test_distTrackFasta_2_out.wig",
		ChromSizes:   "testdata/test_distTrack_1.chrom.sizes",
		ExpectedFile: "testdata/test_distTrackFasta_2_expected.wig",
		DefaultValue: -1,
	},
	
}

func TestDistTrack(t *testing.T) {
	var err error
	var alpha, beta map[string]wig.Wig
	var s DistTrackSettings
	for _, testCase := range distTrackTests {
		s = DistTrackSettings{
			InFile1:       testCase.InFile1,
            InFile2:       testCase.InFile2,
			OutFile:      testCase.OutFile,
			DefaultValue: testCase.DefaultValue,
            OutName: testCase.OutName,
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
