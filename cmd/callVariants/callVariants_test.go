package main

import (
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"os"
	"testing"
)

// TODO add better testing once sam simulation code is done
func TestCallVariants(t *testing.T) {
	var err error
	ref := "testdata/human_chrM.fasta"
	exp := []string{"testdata/human_chrM.bam"}
	norm := []string{"testdata/human_chrM2.bam"}
	outFile := "testdata/test_output.vcf"
	expectedFile := "testdata/test_expected.vcf"
	maxP := 1.1
	minAf := 0.0
	minCoverage := 0
	callVariants(exp, norm, ref, outFile, maxP, minAf, minCoverage)

	if !fileio.AreEqualIgnoreComments(outFile, expectedFile) {
		t.Error("problem with variant calling")
	}

	if !t.Failed() {
		err = os.Remove(outFile)
		exception.PanicOnErr(err)
	}
}
