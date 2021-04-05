package vcf

import (
	"os"
	"testing"
)

var readWriteTests = []struct {
	filename string // input
}{
	{"testdata/test.vcf"},
	{"testdata/GP_DP_Samples.vcf"},
	{"testdata/SingleHapData.vcf"},
}

func TestReadToChan(t *testing.T) {
	alpha := Read("testdata/test.vcf")
	var beta []*Vcf
	vcfPipe, _ := GoReadToChan("testdata/test.vcf")

	for vcfs := range vcfPipe {
		beta = append(beta, vcfs)
	}
	if !AllEqual(alpha, beta) {
		t.Errorf("Error there might be a bug in the ReadToChan implementation\n")
	}
}

func TestWriteAndRead(t *testing.T) {
	var actual []*Vcf
	for _, test := range readWriteTests {
		tempFile := "tmp"
		actual = Read(test.filename)
		Write(tempFile, actual)
		alpha := Read(tempFile)
		beta := Read(test.filename)
		if !AllEqual(alpha, beta) {
			t.Errorf("Error: Read and write files do not match\n")
		}
		os.Remove(tempFile)
	}
}
