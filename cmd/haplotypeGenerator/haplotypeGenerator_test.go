package main

import (
	"os"
	"testing"

	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
)

var haplotypeGeneratorTests = []struct {
	InputFaFile   string
	InputVcfFile  string
	InputBedFile  string
	OutDir        string
	LineLength    int
	Verbose       int
	OutputFiles   []string
	ExpectedFiles []string
}{
	{InputFaFile: "testdata/test.fa",
		InputVcfFile:  "testdata/test.vcf",
		InputBedFile:  "testdata/test.bed",
		OutDir:        "testdata/outdir",
		LineLength:    50,
		Verbose:       0,
		OutputFiles:   []string{"testdata/outdir/CHR1.10.20.fa", "testdata/outdir/CHR1.35.45.fa"},
		ExpectedFiles: []string{"testdata/outdir/expected.CHR1.10.20.fa", "testdata/outdir/expected.CHR1.35.45.fa"}},
}

func TestHaplotypeGenerator(t *testing.T) {
	var err error
	var s Settings
	for _, v := range haplotypeGeneratorTests {
		s = Settings{
			ReferenceGenomeFile: v.InputFaFile,
			VcfFile:             v.InputVcfFile,
			RegionBedFile:       v.InputBedFile,
			OutDir:              v.OutDir,
			LineLength:          v.LineLength,
			Verbose:             v.Verbose,
		}
		haplotypeGenerator(s)
		for i := range v.OutputFiles {
			if !fileio.AreEqual(v.OutputFiles[i], v.ExpectedFiles[i]) {
				t.Errorf("Error in haplotypeGenerator, output was not as expected.")
			} else {
				err = os.Remove(v.OutputFiles[i])
				exception.PanicOnErr(err)
			}
		}
	}
}
