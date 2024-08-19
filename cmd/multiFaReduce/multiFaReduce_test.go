package main

import (
	"github.com/vertgenlab/gonomics/bed"
	"os"
	"testing"

	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fasta"
)

var MfaReduceTests = []struct {
	inputFile    string
	outputFile   string
	expectedFile string
	bedFilename  string
	chrom        string
	expectedBed  string
}{
	{"testdata/test.mfa", "testdata/output.mfa", "testdata/expected.mfa", "", "", ""},
	{"testdata/test.mfa", "testdata/output.mfa", "testdata/expected.mfa", "testdata/output.bed", "chrTest", "testdata/expected.bed"},
	{"testdata/test2.mfa", "testdata/output2.mfa", "testdata/expected2.mfa", "testdata/output2.bed", "chrTest", "testdata/expected2.bed"},
	{"testdata/test3.mfa", "testdata/output3.mfa", "testdata/expected3.mfa", "testdata/output3.bed", "chrTest", "testdata/expected3.bed"},
	{"testdata/test4.mfa", "testdata/output4.mfa", "testdata/expected4.mfa", "testdata.output4.bed", "chrTest", "testdata/expected4.bed"},
}

func TestMfaReduce(t *testing.T) {
	var err error
	for _, v := range MfaReduceTests {
		mfaReduce(v.inputFile, v.outputFile, v.bedFilename, v.chrom)
		records := fasta.Read(v.outputFile)
		expected := fasta.Read(v.expectedFile)
		if !fasta.AllAreEqual(records, expected) {
			t.Errorf("Error in mfaReduce.")
		} else {
			err = os.Remove(v.outputFile)
			exception.PanicOnErr(err)
		}
		if v.bedFilename != "" {
			recordsBed := bed.Read(v.bedFilename)
			expectedBed := bed.Read(v.expectedBed)
			if !bed.AllAreEqual(recordsBed, expectedBed) {
				t.Errorf("Error in mfaReduce bed output.")
			} else {
				err = os.Remove(v.bedFilename)
				exception.PanicOnErr(err)
			}
		}
	}
}
