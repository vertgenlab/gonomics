package main

import (
	"github.com/vertgenlab/gonomics/fileio"
	"os"
	"testing"
)

func TestFaMaskIupac(t *testing.T) {
	in := "testdata/input.fa"
	out := "testdata/test_output.fa"
	expected := "testdata/expected.fa"

	faMaskIupac(in, out)

	if !fileio.AreEqual(out, expected) {
		t.Error("ERROR: problem with faMaskIupac")
	} else {
		os.Remove(out)
	}
}
