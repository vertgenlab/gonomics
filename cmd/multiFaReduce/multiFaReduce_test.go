package main

import (
	"os"
	"testing"

	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fasta"
)

var MfaReduceTests = []struct {
	inputFile    string
	outputFile   string
	expectedFile string
}{
	{"testdata/test.mfa", "testdata/output.mfa", "testdata/expected.mfa"},
}

func TestMfaReduce(t *testing.T) {
	var err error
	for _, v := range MfaReduceTests {
		mfaReduce(v.inputFile, v.outputFile)
		records := fasta.Read(v.outputFile)
		expected := fasta.Read(v.expectedFile)
		if !fasta.AllAreEqual(records, expected) {
			t.Errorf("Error in mfaReduce.")
		} else {
			err = os.Remove(v.outputFile)
			exception.PanicOnErr(err)
		}
	}
}
