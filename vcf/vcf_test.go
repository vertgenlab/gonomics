package vcf

import (
	"log"
	"os"
	"testing"
)

var readWriteTests = []struct {
	filename string // input
}{
	{"testdata/test.vcf"},
}

func TestReadToChan(t *testing.T) {
	alpha := Read("testdata/m54089_181119_141900.subreads.vcf")
	var beta []*Vcf
	vcfPipe := make(chan *Vcf)
	go ReadToChan("testdata/m54089_181119_141900.subreads.vcf", vcfPipe)
	for vcfs := range vcfPipe {
		beta = append(beta, vcfs)
	}
	log.Printf("alpha=%d, beta=%d", len(alpha), len(beta))
	if len(alpha) != len(beta) {
		t.Errorf("Vcf lengths are not equal...")
	}
	if !AllEqual(alpha, beta) {
		t.Errorf("VCF files are not the same...")
	}
	PrintVcfLines(alpha, 5)
}

func TestWriteAndRead(t *testing.T) {
	var actual []*Vcf
	for _, test := range readWriteTests {
		tempFile := "tmp"
		actual = Read(test.filename)
		Write(tempFile, actual)
		alpha := Read(tempFile)
		beta := Read("testdata/test.vcf")
		PrintVcfLines(beta, 5)
		log.Printf("alpha=%d, beta=%d", len(alpha), len(beta))
		if len(alpha) != len(beta) {
			t.Errorf("Vcf lengths are not equal...")
		}
		if !AllEqual(alpha, beta) {
			t.Errorf("VCF files are not the same...")
		}
		os.Remove(tempFile)
	}
}
