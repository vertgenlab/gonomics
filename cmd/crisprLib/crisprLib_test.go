package main

import (
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"os"
	"testing"
)

func TestCrisprLib(t *testing.T) {
	inFa := "testdata/inGuides.fa"
	exp := "testdata/exp.oligo.fa"
	out := "testdata/out.oligo.fa"

	crisprLib(inFa, out, 2)

	if !fileio.AreEqual(exp, out) {
		t.Errorf("ERROR in crisprLib. Expected file isn't equal to the output.")
	} else {
		err := os.Remove(out)
		exception.PanicOnErr(err)
	}
}
