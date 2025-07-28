package main

import (
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"os"
	"testing"
)

var FastqReplaceSeqTests = []struct {
	inFile              string
	outFile             string
	expectedFile        string
	findReplaceFile     string
	findReplaceDelim    string
	ignoreCase          bool
	replacedRecordsOnly bool
}{
	{"testdata/test1.fastq", "tmpOut.fastq", "testdata/expected_allFile.fastq", "testdata/findReplace1.txt", "\t", false, false},
	{"testdata/test1.fastq", "tmpOut.fastq", "testdata/expected_onlyReplaced.fastq", "testdata/findReplace1.txt", "\t", false, true},
}

func TestFastqReplaceSeq(t *testing.T) {
	var err error
	var s Settings

	for _, v := range FastqReplaceSeqTests {
		s = Settings{
			inFile:              v.inFile,
			outFile:             v.outFile,
			findReplaceFile:     v.findReplaceFile,
			findReplaceDelim:    v.findReplaceDelim,
			ignoreCase:          v.ignoreCase,
			replacedRecordsOnly: v.replacedRecordsOnly,
		}
		fastqReplaceSeq(s)
		if !fileio.AreEqual(v.outFile, v.expectedFile) {
			t.Errorf("Error in findAndReplace")
		} else {
			err = os.Remove(v.outFile)
			exception.PanicOnErr(err)
		}
	}

}
