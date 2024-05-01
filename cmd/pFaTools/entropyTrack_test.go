package main

import (
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/wig"
	"os"
	"testing"
)

var entropyTrackTests = []struct {
	InFile       string
	OutFile      string
	ChromSizes   string
	ExpectedFile string
	DefaultValue float64
}{
	{InFile: "testdata/test_sample_input.pfa",
		OutFile:      "testdata/entropyTrack.wig",
		ChromSizes:   "testdata/test_sample_input.chrom.sizes",
		ExpectedFile: "testdata/expected.EntropyTrack.wig",
		DefaultValue: -1,
	},
}

func TestEntropyTrack(t *testing.T) {
	var err error
	var alpha, beta map[string]wig.Wig
	var s EntropyTrackSettings
	for _, testCase := range entropyTrackTests {
		s = EntropyTrackSettings{
			InFile:       testCase.InFile,
			OutFile:      testCase.OutFile,
			DefaultValue: testCase.DefaultValue,
		}
		pFaEntropyTrack(s)
		alpha = wig.Read(testCase.OutFile, testCase.ChromSizes, testCase.DefaultValue)
		beta = wig.Read(testCase.ExpectedFile, testCase.ChromSizes, testCase.DefaultValue)
		if !wig.AllEqual(alpha, beta, 1e-6) {
			t.Errorf("Error: EntropyTrack output was not as expected.\n")
		} else {
			err = os.Remove(testCase.OutFile)
			exception.PanicOnErr(err)
		}
	}
}
