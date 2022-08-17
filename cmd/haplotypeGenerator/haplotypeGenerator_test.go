package main

import (
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"os"
	"testing"
)

var haplotypeGeneratorTests = []struct {
	inputFaFile   string
	inputVcfFile  string
	inputBedFile  string
	outdir        string
	outputFiles   []string
	expectedFiles []string
}{
	{"testdata/test.fa", "testdata/test.vcf", "testdata/test.bed", "testdata/outdir", []string{"testdata/outdir/CHR1.10.20.fa", "testdata/outdir/CHR1.35.45.fa"}, []string{"testdata/outdir/expected.CHR1.10.20.fa", "testdata/outdir/expected.CHR1.35.45.fa"}},
}

func TestHaplotypeGenerator(t *testing.T) {
	var err error
	for _, v := range haplotypeGeneratorTests {
		haplotypeGenerator(v.inputFaFile, v.inputVcfFile, v.inputBedFile, v.outdir)
		for i := range v.outputFiles {
			if !fileio.AreEqual(v.outputFiles[i], v.expectedFiles[i]) {
				t.Errorf("Error in haplotypeGenerator, output was not as expected.")
			} else {
				for i := range v.outputFiles {
					err = os.Remove(v.outputFiles[i])
					exception.PanicOnErr(err)
				}
			}
		}

	}
}
