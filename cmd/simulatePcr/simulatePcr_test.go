package main

import (
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"os"
	"testing"
)

func TestSimPcr(t *testing.T) {
	template := "testdata/test.fasta"
	primers := []string{
		"GCCTCCGTGAGGCTAC",
		"TTGAGGATCTTTTCTTCACG",
	}
	simulatePcr(primers, template, "testdata/actual1.bed", 1000)

	if !fileio.AreEqual("testdata/actual1.bed", "testdata/expected1.bed") {
		t.Error("problem with simulatePcr")
	}

	if !t.Failed() {
		err := os.Remove("testdata/actual1.bed")
		exception.PanicOnErr(err)
	}

	primers = []string{
		"ATG",
	}
	simulatePcr(primers, template, "testdata/actual2.bed", 1000)

	if !fileio.AreEqual("testdata/actual2.bed", "testdata/expected2.bed") {
		t.Error("problem with simulatePcr")
	}

	if !t.Failed() {
		err := os.Remove("testdata/actual2.bed")
		exception.PanicOnErr(err)
	}
}
