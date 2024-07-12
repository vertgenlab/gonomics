package main

import (
	"testing"

	"github.com/vertgenlab/gonomics/fileio"
)

var pcrSimTests = []struct {
	Primers       []string
	BedFile       string
	ExpectedBed   string
	FastqFile     string
	ExpectedFastq string
	Length        int
	WithPrimer    bool
}{
	{
		Primers: []string{
			"GCCTCCGTGAGGCTAC",
			"TTGAGGATCTTTTCTTCACG",
		},
		BedFile:       "testdata/actual1.bed",
		ExpectedBed:   "testdata/expected1.bed",
		FastqFile:     "testdata/actual1.fastq",
		ExpectedFastq: "testdata/expected1.fastq",
		Length:        1000,
		WithPrimer:    false,
	},
	{Primers: []string{
		"ATG",
	},
		BedFile:       "testdata/actual2.bed",
		ExpectedBed:   "testdata/expected2.bed",
		FastqFile:     "testdata/actual2.fastq",
		ExpectedFastq: "testdata/expected2.fastq",
		Length:        1000,
		WithPrimer:    true,
	},
}

func TestSimPcr(t *testing.T) {
	template := "testdata/test.fasta"
	for _, pcr := range pcrSimTests {
		simulatePcr(pcr.Primers, template, pcr.BedFile, pcr.FastqFile, pcr.Length, pcr.WithPrimer)
		if !fileio.AreEqual(pcr.ExpectedBed, pcr.ExpectedFastq) {
			t.Error("Error: Problem with simulatePcr")
		}
		if !t.Failed() {
			fileio.EasyRemove(pcr.BedFile)
			fileio.EasyRemove(pcr.FastqFile)
		}
	}
}
