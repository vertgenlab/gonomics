package main

import (
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"os"
	"testing"
)

func TestFaUniq(t *testing.T) {
	inputFile := "testdata/test_in.fa"
	outputFile := "testdata/test_out.fa"
	expectedFile := "testdata/expected_out.fa"

	var err error
	faUniq(inputFile, outputFile)
	if !fileio.AreEqual(outputFile, expectedFile) {
		t.Errorf("Error in faUniq, output was not as expected.")
	} else {
		err = os.Remove(outputFile)
		exception.PanicOnErr(err)
	}
}