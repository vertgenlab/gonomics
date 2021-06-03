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
	alpha, _ := Read("testdata/test.vcf")
	var beta []Vcf
	vcfPipe, _ := GoReadToChan("testdata/test.vcf")

	for vcfs := range vcfPipe {
		beta = append(beta, vcfs)
	}
	if !AllEqual(alpha, beta) {
		t.Errorf("Error there might be a bug in the ReadToChan implementation\n")
	}
}

func TestWriteAndRead(t *testing.T) {
	var actual []Vcf
	for _, test := range readWriteTests {
		tempFile := "tmp"
		actual, _ = Read(test.filename)
		Write(tempFile, actual)
		alpha, _ := Read(tempFile)
		beta, _ := Read(test.filename)
		if !AllEqual(alpha, beta) {
			t.Errorf("Error: Read and write files do not match\n")
		}
		os.Remove(tempFile)
	}
}

func BenchmarkRead(b *testing.B) {
	for i := 0; i < b.N; i++ {
		data, _ := GoReadToChan("testdata/test.vcf")
		for _ = range data { // stall
		}
	}
}

func BenchmarkWrite(b *testing.B) {
	records, _ := Read("testdata/test.vcf")
	for i := 0; i < b.N; i++ {
		Write("tmp", records)
	}
	os.Remove("tmp")
}
