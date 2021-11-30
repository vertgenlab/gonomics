package main

import (
	"testing"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"os"
)

var VcfAncestorAnnotationTests = []struct {
	Infile string
	Outfile string
	ExpectedFile string
	FaFile string
}{
	{"testdata/in.vcf", "testdata/tmpOut.vcf", "testdata/expected.vcf", "testdata/test.fa"},
}

func TestVcfAncestorAnnotation(t *testing.T) {
	var err error
	for _, v := range VcfAncestorAnnotationTests {
		vcfAncestorAnnotation(v.Infile, v.FaFile, v.Outfile)
		if !fileio.AreEqual(v.Outfile, v.ExpectedFile) {
			t.Errorf("Error in VcfAncestorAnnotation.")
		} else {
			err = os.Remove(v.Outfile)
			exception.PanicOnErr(err)
		}
	}
}