package main

import (
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"os"
	"sort"
	"testing"
)

var QuickOrthologTests = []struct {
	TName                string
	QName                string
	GtfFile              string
	ChainFile            string
	ChromSizes           string
	OutFile              string
	Unmapped             string
	CanonicalTranscript  bool
	ExpectedOutFile      string
	ExpectedUnmappedFile string
}{
	{TName: "hg38",
		QName:                "panTro6",
		GtfFile:              "testdata/chrM.hg38.panTro6.gtf",
		ChainFile:            "testdata/chrM.hg38.panTro6.chain",
		ChromSizes:           "testdata/chrM.chrom.sizes",
		OutFile:              "testdata/tmp.Out.txt",
		Unmapped:             "testdata/tmp.unmapped.txt",
		CanonicalTranscript:  false,
		ExpectedOutFile:      "testdata/expected.out.txt",
		ExpectedUnmappedFile: "testdata/expected.unmapped.txt",
	},
}

func TestQuickOrthologs(t *testing.T) {
	var err error
	var records, expectedRecords []string
	var s Settings
	for _, v := range QuickOrthologTests {
		s = Settings{
			TName:               v.TName,
			QName:               v.QName,
			GtfFile:             v.GtfFile,
			ChainFile:           v.ChainFile,
			ChromSizes:          v.ChromSizes,
			OutFile:             v.OutFile,
			Unmapped:            v.Unmapped,
			CanonicalTranscript: v.CanonicalTranscript,
		}
		quickOrthologs(s)
		//we have to sort the output because of random hash order
		records = fileio.Read(v.OutFile)
		sort.Strings(records)
		expectedRecords = fileio.Read(v.ExpectedOutFile)
		sort.Strings(expectedRecords)

		if !allAreEqual(records, expectedRecords) {
			t.Errorf("Error in quickOrthologs.")
		} else {
			err = os.Remove(v.OutFile)
			exception.PanicOnErr(err)
		}

		records = fileio.Read(v.Unmapped)
		sort.Strings(records)
		expectedRecords = fileio.Read(v.ExpectedUnmappedFile)
		sort.Strings(expectedRecords)

		if !allAreEqual(records, expectedRecords) {
			t.Errorf("Error in quickOrthologs.")
		} else {
			err = os.Remove(v.Unmapped)
			exception.PanicOnErr(err)
		}
	}
}

func allAreEqual(a []string, b []string) bool {
	if len(a) != len(b) {
		return false
	}
	for i := range a {
		if a[i] != b[i] {
			return false
		}
	}
	return true
}
