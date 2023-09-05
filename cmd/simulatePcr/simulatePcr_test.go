package main

import (
	"os"
	"testing"

	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
)

func TestSimPcr(t *testing.T) {
	template := "testdata/test.fasta"
	primers := []string{
		"GCCTCCGTGAGGCTAC",
		"TTGAGGATCTTTTCTTCACG",
	}
	simulatePcr(primers, template, "testdata/actual1.bed", "testdata/actual1.fastq", 1000, false)

	if !fileio.AreEqual("testdata/actual1.bed", "testdata/expected1.bed") || !fileio.AreEqual("testdata/actual1.fastq", "testdata/expected1.fastq") {
		t.Error("problem with simulatePcr")
	}

	var err error
	if !t.Failed() {
		err = os.Remove("testdata/actual1.bed")
		exception.PanicOnErr(err)
		err = os.Remove("testdata/actual1.fastq")
		exception.PanicOnErr(err)
	}

	primers = []string{
		"ATG",
	}
	simulatePcr(primers, template, "testdata/actual2.bed", "testdata/actual2.fastq", 1000, true)

	if !fileio.AreEqual("testdata/actual2.bed", "testdata/expected2.bed") || !fileio.AreEqual("testdata/actual2.fastq", "testdata/expected2.fastq") {
		t.Error("problem with simulatePcr")
	}

	if !t.Failed() {
		err = os.Remove("testdata/actual2.bed")
		exception.PanicOnErr(err)
		err = os.Remove("testdata/actual2.fastq")
		exception.PanicOnErr(err)
	}
}
