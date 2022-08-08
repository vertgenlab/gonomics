package main

import (
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"os"
	"testing"
)

var divergenceSpectrumTests = []struct {
	InFile       string
	VariantsFile string
	OutFile      string
	ExpectedFile string
}{
	{"testdata/test.bed", "testdata/test.vcf", "testdata/tmp.bed", "testdata/expected.bed"},
}

func TestDivergenceSpectrum(t *testing.T) {
	var err error
	for _, v := range divergenceSpectrumTests {
		divergenceSpectrum(v.InFile, v.VariantsFile, v.OutFile)
		if !fileio.AreEqual(v.OutFile, v.ExpectedFile) {
			t.Errorf("Error in divergenceSpectrum. Output not as expected.")
		} else {
			err = os.Remove(v.OutFile)
			exception.PanicOnErr(err)
		}
	}
}
