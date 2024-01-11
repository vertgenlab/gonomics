package main

import (
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"os"
	"testing"
)

var SamToWigTests = []struct {
	InFile       string
	Reference    string
	OutFile      string
	ExpectedFile string
	FragLength   int
	DefaultValue float64
	Deletions    bool
}{
	{InFile: "testdata/test1.sam",
		Reference:    "testdata/test.chrom.sizes",
		OutFile:      "testdata/tmp.test1.wig",
		ExpectedFile: "testdata/test1.wig",
		FragLength:   -1,
		DefaultValue: 0,
		Deletions:    false},
	{InFile: "testdata/test2.sam",
		Reference:    "testdata/test.chrom.sizes",
		OutFile:      "testdata/tmp.test2.wig",
		ExpectedFile: "testdata/test2.wig",
		FragLength:   30,
		DefaultValue: 0,
		Deletions:    false},
	{InFile: "testdata/test1.bam",
		Reference:    "testdata/test.chrom.sizes",
		OutFile:      "testdata/tmp.bam.wig",
		ExpectedFile: "testdata/test1.wig",
		FragLength:   -1,
		DefaultValue: 0,
		Deletions:    false},
	{InFile: "testdata/test2.bam",
		Reference:    "testdata/test.chrom.sizes",
		OutFile:      "testdata/tmp.bamFrag.wig",
		ExpectedFile: "testdata/test2.wig",
		FragLength:   30,
		DefaultValue: 0,
		Deletions:    false},
	{InFile: "testdata/test1.sam",
		Reference:    "testdata/test.chrom.sizes",
		OutFile:      "testdata/tmp.withDel.wig",
		ExpectedFile: "testdata/test1.withDel.wig",
		FragLength:   -1,
		DefaultValue: 0,
		Deletions:    true},
}

func TestSamToWig(t *testing.T) {
	var s Settings
	for _, v := range SamToWigTests {
		s = Settings{
			SamFileName:        v.InFile,
			ChromSizesFileName: v.Reference,
			OutFileName:        v.OutFile,
			FragLength:         v.FragLength,
			DefaultValue:       v.DefaultValue,
			Deletions:          v.Deletions,
		}
		samToWig(s)
		if !fileio.AreEqual(v.OutFile, v.ExpectedFile) {
			t.Errorf("Error in samToWig")
		} else {
			err := os.Remove(v.OutFile)
			exception.PanicOnErr(err)
		}
	}
}
