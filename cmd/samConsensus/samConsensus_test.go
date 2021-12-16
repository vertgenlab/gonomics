package main

import (
	"github.com/vertgenlab/gonomics/common"
	"github.com/vertgenlab/gonomics/fileio"
	"os"
	"testing"
)

var samConsensusTests = []struct {
	inFile           string
	refFile          string
	outFile_expected string
	vcfFile_expected string
}{
	{"testdata/test.sam", "testdata/test.ref.fa", "testdata/test.out.fa", "testdata/test.out.vcf"},
}

func TestsamConsensus(t *testing.T) {
	for _, v := range samConsensusTests {
		samConsensus(v.inFile, v.refFile, "outFile_tmp.fa", "vcfFile_tmp.vcf")
		if !fileio.AreEqual("outFile_tmp.fa", v.outFile_expected) {
			t.Errorf("Error in samConsensus: generating output fa file")
		}
		err := os.Remove("outFile_tmp.fa")
		if err != nil {
			common.ExitIfError(err)
		}
		if !fileio.AreEqual("vcfFile_tmp.vcf", v.vcfFile_expected) {
			t.Errorf("Error in samConsensus: generating output vcf file")
		}
		err = os.Remove("vcfFile_tmp.vcf")
		if err != nil {
			common.ExitIfError(err)
		}
	}
}
