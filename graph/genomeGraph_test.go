package graph

import (
	"github.com/vertgenlab/gonomics/axt"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/vcf"
	"os"
	"testing"
)

var readWriteTests = []struct {
	filename string // input
}{
	{"testdata/dev.gg"},
}

func TestRead(t *testing.T) {
	for _, test := range readWriteTests {
		_ = Read(test.filename)
	}
}

func TestAxtToGraph(t *testing.T) {
	axtAlign := axt.Read("testdata/test.axt")
	vcfAxt := make([]*vcf.Vcf, 0)
	for i := 0; i < len(axtAlign); i++ {
		vcfAxt = append(vcfAxt, axt.AxtToVcf(axtAlign[i])...)
	}
	ref := fasta.Read("testdata/test.fa")
	g := RefernceToGraph(vcfAxt, ref, NewGraph())
	Write("temp.gg", g)
	gsw := Read("temp.gg")
	PrintGraph(gsw)
	err := os.Remove("temp.gg")
	if err != nil {
		t.Errorf("Deleting temp file %s gave an error.", "temp.gg")
	}
}

func TestWriteAndRead(t *testing.T) {
	var actual *GenomeGraph
	for _, test := range readWriteTests {
		tempFile := test.filename + ".tmp"
		actual = Read(test.filename)
		Write(tempFile, actual)
		err := os.Remove(tempFile)
		if err != nil {
			t.Errorf("Deleting temp file %s gave an error.", tempFile)
		}
	}
}
