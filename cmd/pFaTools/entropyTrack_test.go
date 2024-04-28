package main

import (
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"os"
	"testing"
)

var entropyTrackTests = []struct {
	InFile       string
	OutFile      string
	ExpectedFile string
	DefaultValue float64
}{
	{InFile: "testdata/test_sample_input.pfa",
		OutFile:      "testdata/entropyTrack.wig",
		ExpectedFile: "testdata/expected.EntropyTrack.wig",
		DefaultValue: -1,
	},
}

func TestEntropyTrack(t *testing.T) {
	var err error
	var s EntropyTrackSettings
	for _, testCase := range entropyTrackTests {
		s = EntropyTrackSettings{
			InFile:       testCase.InFile,
			OutFile:      testCase.OutFile,
			DefaultValue: testCase.DefaultValue,
		}
		pFaEntropyTrack(s)
		if !fileio.AreEqual(testCase.OutFile, testCase.ExpectedFile) {
			t.Errorf("Error: EntropyTrack output was not as expected.\n")
		} else {
			err = os.Remove(testCase.OutFile)
			exception.PanicOnErr(err)
		}
	}
}
