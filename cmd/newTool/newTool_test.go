package main

import (
	"github.com/vertgenlab/gonomics/fileio"
	"log"
	"os"
	"strings"
	"testing"
)

func TestNewTool(t *testing.T) {
	functionForTheTool("testdata/test.bed", "testdata/test.fasta", "testdata/test.out")
	testOut := fileio.Read("testdata/test.out")
	expectedOut := fileio.Read("testdata/expected.out")

	for i := range testOut {
		if strings.Compare(testOut[i], expectedOut[i]) != 0 {
			log.Fatalf("Expected and test did not match. Expected: %s, Test: %s", expectedOut, testOut)
		}
	}
	os.Remove("testdata/test.out")
}
