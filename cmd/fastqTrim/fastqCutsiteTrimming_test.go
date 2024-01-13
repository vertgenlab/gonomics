package main

import (
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"os"
	"testing"
)

func TestCutSiteTrim(t *testing.T) {
	in := "testdata/in.fastq"
	bases := []dna.Base{2, 0, 3, 1}
	out := "testdata/out.fastq"
	exp := "testdata/exp.fastq"

	cutSiteTrim(in, out, bases, 0)

	if !fileio.AreEqual(out, exp) {
		t.Errorf("Error in CutSiteTrim")
	} else {
		err := os.Remove(out)
		exception.PanicOnErr(err)
	}
}
