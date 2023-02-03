package main

import (
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"math/rand"
	"os"
	"testing"
)

var samConsensusTests = []struct {
	inFile             string
	refFile            string
	outFile            string
	vcfFile            string
	chainFile          string
	substitutionsOnly  bool
	insertionThreshold float64
	outFile_expected   string
	vcfFile_expected   string
}{
	{inFile: "testdata/test.sam",
		refFile:            "testdata/test.ref.fa",
		outFile:            "testdata/tmpOutFile.fa",
		vcfFile:            "testdata/tmpVcfFile.vcf",
		chainFile:          "",
		substitutionsOnly:  true,
		insertionThreshold: 0.9,
		outFile_expected:   "testdata/test.out.fa",
		vcfFile_expected:   "testdata/test.out.vcf"},
	{inFile: "testdata/test.sam",
		refFile:            "testdata/test.ref.fa",
		outFile:            "testdata/tmpOutFile.indel.fa",
		vcfFile:            "testdata/tmpVcfFile.indel.vcf",
		chainFile:          "",
		substitutionsOnly:  false,
		insertionThreshold: 0.9,
		outFile_expected:   "testdata/test.out.indel.fa",
		vcfFile_expected:   "testdata/test.out.indel.vcf"},
}

func TestSamConsensus(t *testing.T) {
	rand.Seed(1)
	var err error
	var s Settings
	for _, v := range samConsensusTests {
		s = Settings{
			SamFileName:        v.inFile,
			RefFile:            v.refFile,
			OutFile:            v.outFile,
			VcfFile:            v.vcfFile,
			ChainFile:          v.chainFile,
			SubstitutionsOnly:  v.substitutionsOnly,
			InsertionThreshold: v.insertionThreshold,
		}
		samConsensus(s)

		if !fileio.AreEqual(v.outFile, v.outFile_expected) {
			t.Errorf("Error in samConsensus: generating output fa file")
		} else {
			err = os.Remove(v.outFile)
			exception.PanicOnErr(err)
		}

		if v.vcfFile != "" && !fileio.AreEqual(v.vcfFile, v.vcfFile_expected) {
			t.Errorf("Error in samConsensus: generating output vcf file")
		} else {
			err = os.Remove(v.vcfFile)
			exception.PanicOnErr(err)
		}
	}
}
