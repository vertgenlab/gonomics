package main

import (
	"github.com/vertgenlab/gonomics/fileio"
	"log"
	"os"
	"strings"
	"testing"
)

func TestExampleTool(t *testing.T) { //compare the output results made from the predetermined files and see if the program gives th expected result
	bedToAminoAcid("testdata/test.bed", "testdata/test.fasta", "testdata/test.out")
	testOut := fileio.Read("testdata/test.out")
	expectedOut := fileio.Read("testdata/expected.txt")

	for i := range testOut { //check if the output is equal to the expected result
		if strings.Compare(testOut[i], expectedOut[i]) != 0 {
			log.Fatalf("Expected and test did not match. Expected: %s, Test: %s", expectedOut, testOut) //if the files aren't equal write an error
		}
	}
	os.Remove("testdata/test.out")
}
