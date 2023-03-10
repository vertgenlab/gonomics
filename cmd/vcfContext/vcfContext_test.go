package main

import (
	"github.com/vertgenlab/gonomics/fileio"
	"os"
	"testing"
)

func TestVcfContext(t *testing.T) {
	fa := "testdata/test.fasta"
	v := "testdata/test.vcf"
	expectedMerge := "testdata/expectedMergeComplements.txt"
	expectedInclude := "testdata/expectedIncludeComplements.txt"
	outMerge := "testdata/outMergeComplements.txt"
	outInclude := "testdata/outIncludeComplements.txt"

	vcfContext(fa, v, outMerge, 1, false, 0)
	vcfContext(fa, v, outInclude, 1, true, 0)

	if !fileio.AreEqual(expectedMerge, outMerge) || !fileio.AreEqual(expectedInclude, outInclude) {
		t.Error("ERROR: problem with vcfContext")
		return
	}

	os.Remove(outMerge)
	os.Remove(outInclude)
}
