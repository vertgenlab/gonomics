package main

import (
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"os"
	"testing"
)

func TestBedTrim(t *testing.T) {
	bedTrim(30, "testdata/in.bed", "testdata/out.30.bed")
	if !fileio.AreEqual("testdata/exp.30.bed", "testdata/out.30.bed") {
		t.Errorf("Error in bedTrim: expected and observed output files are not equal.")
	} else {
		exception.PanicOnErr(os.Remove("testdata/out.30.bed"))
	}
	bedTrim(100, "testdata/in.bed", "testdata/out.100.bed")
	if !fileio.AreEqual("testdata/exp.100.bed", "testdata/out.100.bed") {
		t.Errorf("Error in bedTrim: expected and observed output files are not equal.")
	} else {
		exception.PanicOnErr(os.Remove("testdata/out.100.bed"))
	}
}
