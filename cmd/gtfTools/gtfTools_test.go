package main

import (
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"os"
	"testing"
)

var ToBedTests = []struct {
	InFile         string
	OutFile        string
	ExpectedFile   string
	Tss            bool
	ChromSizesFile string
	Merge          bool
}{
	{InFile: "testdata/test.gtf",
		OutFile:        "testdata/tmp.bed",
		ExpectedFile:   "testdata/testOut.bed",
		Tss:            false,
		ChromSizesFile: "",
		Merge:          false,
	},
	{InFile: "testdata/test.gtf",
		OutFile:        "testdata/tmp.tss.bed",
		ExpectedFile:   "testdata/expected.tss.bed",
		Tss:            true,
		ChromSizesFile: "testdata/chr1.chrom.sizes",
		Merge:          false,
	},
}

func TestToBed(t *testing.T) {
	var err error
	var s ToBedSettings
	for _, v := range ToBedTests {
		s = ToBedSettings{
			InFile:        v.InFile,
			OutFile:       v.OutFile,
			Tss:           v.Tss,
			ChromSizeFile: v.ChromSizesFile,
			Merge:         v.Merge,
		}
		toBed(s)
		if !fileio.AreEqual(v.OutFile, v.ExpectedFile) {
			t.Errorf("Error: in gtfTools toBed, output was not as expected.")
		} else {
			err = os.Remove(v.OutFile)
			exception.PanicOnErr(err)
		}
	}
}

var FilterTests = []struct {
	InFile           string
	OutFile          string
	ExpectedFile     string
	GeneNameList     string
	ChromFilter      string
	CodingTranscript bool
}{
	{InFile: "../../gtf/testdata/test.gtf",
		OutFile:          "testdata/tmp.filter.gtf",
		ExpectedFile:     "testdata/expected.filter.gtf",
		GeneNameList:     "testdata/geneList.txt",
		ChromFilter:      "",
		CodingTranscript: false,
	},
	{InFile: "testdata/chromFilter.gtf",
		OutFile:          "testdata/tmp.filter.gtf",
		ExpectedFile:     "testdata/expected.chromFilter.gtf",
		GeneNameList:     "",
		ChromFilter:      "chrM",
		CodingTranscript: false,
	},
	{InFile: "testdata/chromFilter.gtf",
		OutFile:          "testdata/tmp.filter.gtf",
		ExpectedFile:     "testdata/expected.chromFilterGeneFilter.gtf",
		GeneNameList:     "testdata/geneListForChromFilter.txt",
		ChromFilter:      "chr1",
		CodingTranscript: false,
	},
	{InFile: "testdata/codingTranscriptFilter.gtf",
		OutFile:          "testdata/tmp.filter.gtf",
		ExpectedFile:     "testdata/expected.codingTranscriptFilter.gtf",
		GeneNameList:     "",
		ChromFilter:      "",
		CodingTranscript: true,
	},
}

func TestFilter(t *testing.T) {
	var s FilterSettings
	var err error
	for _, v := range FilterTests {
		s = FilterSettings{
			InFile:           v.InFile,
			OutFile:          v.OutFile,
			GeneNameList:     v.GeneNameList,
			ChromFilter:      v.ChromFilter,
			CodingTranscript: v.CodingTranscript,
		}
		gtfFilter(s)
		if !fileio.AreEqual(v.OutFile, v.ExpectedFile) {
			t.Errorf("Error: gtfTools filter - outfile was not as expected.")
		} else {
			err = os.Remove(v.OutFile)
			exception.PanicOnErr(err)
		}
	}
}
