package main

import (
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"os"
	"testing"
)

inputFile := test_in.fa
outputFile := test_out.fa
expectedFile := expected_out.fa

func TestFaUniq(t *testing.T) {
	var err error
	faUniq(inputFile, outputFile)
	if !fileio.AreEqual(outputFile, expectedFile) {
		t.Errorf("Error in faUniq, output was not as expected.")
	} else {
		err = os.Remove(outputFile)
		exception.PanicOnErr(err)
	}
}