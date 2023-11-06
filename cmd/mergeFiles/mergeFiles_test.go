package main

import (
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"os"
	"testing"
)

func TestMergeFiles(t *testing.T) {
	in1 := "testdata/file1.sam"
	in2 := "testdata/file2.sam"
	exp := "testdata/exp.sam"
	out := "testdata/out.sam"

	mergeFiles(in1, in2, out)

	if !fileio.AreEqual(exp, out) {
		t.Errorf("expected and output files are not equal")
	} else {
		err := os.Remove(out)
		exception.PanicOnErr(err)
	}
}
