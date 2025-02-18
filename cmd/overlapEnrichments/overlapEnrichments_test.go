package main

import (
	"os"
	"testing"

	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
)

var OverlapEnrichmentsTests = []struct {
	Method            string
	Elements1File     string
	Elements2File     string
	SearchSpaceFile   string
	OutFile           string
	ExpectedFile      string
	TrimToSearchSpace bool
	Verbose           int
	SecondFileList    string
	Relationship      string
}{
	{Method: "exact",
		Elements1File:     "testdata/elements1.bed",
		Elements2File:     "testdata/elements2.bed",
		SearchSpaceFile:   "testdata/tinyNoGap.bed",
		OutFile:           "testdata/tmp.txt",
		ExpectedFile:      "testdata/elements1.elements2.enrichment.txt",
		TrimToSearchSpace: false,
		Verbose:           0,
		SecondFileList:    "",
		Relationship:      "within",
	},
	{Method: "exact",
		Elements1File:     "testdata/elements1.bed",
		Elements2File:     "testdata/elements1.bed",
		SearchSpaceFile:   "testdata/tinyNoGap.bed",
		ExpectedFile:      "testdata/elements1.elements1.enrichment.txt",
		OutFile:           "testdata/self.txt",
		TrimToSearchSpace: false,
		Verbose:           0,
		SecondFileList:    "",
		Relationship:      "within",
	},
	{Method: "exact",
		Elements1File:     "testdata/elements1.bed",
		Elements2File:     "testdata/elements2.bed",
		SearchSpaceFile:   "testdata/tinyNoGap.bed",
		OutFile:           "testdata/trim.txt",
		ExpectedFile:      "testdata/elements1.elements2.enrichment.txt",
		TrimToSearchSpace: true,
		Verbose:           0,
		SecondFileList:    "",
		Relationship:      "within",
	}, //should have no effect to trim to ref Genome if all elements are in the genome.
	{Method: "exact",
		Elements1File:     "testdata/elements1.bed",
		Elements2File:     "testdata/elements3.bed",
		SearchSpaceFile:   "testdata/tinyNoGap.bed",
		OutFile:           "testdata/trim.outside.txt",
		ExpectedFile:      "testdata/elements1.elements3.enrichment.txt",
		Verbose:           0,
		SecondFileList:    "",
		Relationship:      "within",
		TrimToSearchSpace: true}, //elements3 is elements2 with extra elements outside the genome, should be the same answer as the previous check.
	{Method: "exact",
		Elements1File:     "testdata/elements1.bed",
		Elements2File:     "testdata/elements1.bed", // this will be ignored due to the SecondFileList option.
		SearchSpaceFile:   "testdata/tinyNoGap.bed",
		OutFile:           "testdata/test.listOfFiles.enrichment.txt",
		ExpectedFile:      "testdata/expected.listOfFiles.txt",
		Verbose:           0,
		SecondFileList:    "testdata/listOfFiles.txt",
		TrimToSearchSpace: true,
		Relationship:      "within",
	},
	{Method: "exact",
		Elements1File:     "testdata/elements1.bed",
		Elements2File:     "testdata/elements3.bed",
		SearchSpaceFile:   "testdata/tinyNoGap.bed",
		OutFile:           "testdata/trim.outside.any.txt",
		ExpectedFile:      "testdata/elements1.elements3.enrichment.any.txt",
		Verbose:           0,
		SecondFileList:    "",
		Relationship:      "any",
		TrimToSearchSpace: true},
}

func TestOverlapEnrichments(t *testing.T) {
	var err error
	var s Settings
	for _, v := range OverlapEnrichmentsTests {
		s = Settings{
			Method:            v.Method,
			InFile:            v.Elements1File,
			SecondFile:        v.Elements2File,
			SearchSpaceFile:   v.SearchSpaceFile,
			OutFile:           v.OutFile,
			Verbose:           v.Verbose,
			TrimToSearchSpace: v.TrimToSearchSpace,
			SecondFileList:    v.SecondFileList,
			Relationship:      v.Relationship,
		}
		overlapEnrichments(s)
		if !fileio.AreEqual(v.OutFile, v.ExpectedFile) {
			t.Errorf("Error in overlapEnrichments.")
		} else {
			err = os.Remove(v.OutFile)
			exception.PanicOnErr(err)
		}
	}
}
