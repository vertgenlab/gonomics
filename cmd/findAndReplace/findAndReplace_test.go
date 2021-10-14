package main

import (
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"os"
	"testing"
)

var findAndReplaceTestsColumnSpecific = []struct {
	inFile				string
	findAndReplaceFile	string
	outFile				string
	expectedFile		string
	columnNumber		int

}{
	{"testdata/inputFileFake.tsv", "testdata/findReplaceFake.tsv", "testdata/outputFileCreatedColumn0.tsv",
		"testdata/outputFileExpectedColumn0.tsv", 0},
	{"testdata/inputFileFake.tsv", "testdata/findReplaceFake.tsv", "testdata/outputFileCreatedColumn1.tsv",
		"testdata/outputFileExpectedColumn1.tsv", 1},
	{"testdata/inputFileFake.tsv", "testdata/findReplaceFake.tsv", "testdata/outputFileCreatedColumn2.tsv",
		"testdata/outputFileExpectedColumn2.tsv", 2},

}
func TestFindAndReplaceColumnSpecific(t *testing.T) {
	var err error
	for _, v := range findAndReplaceTestsColumnSpecific{
		findAndReplace(v.inFile, v.findAndReplaceFile, v.outFile, v.columnNumber)
		if !fileio.AreEqual(v.outFile, v.expectedFile){
			t.Errorf("Error in findAndReplace")
		}else{
			err = os.Remove(v.outFile)
			exception.PanicOnErr(err)
		}
	}
}



var findAndReplaceTestWholeFileScan = []struct {
	inFile				string
	findAndReplaceFile	string
	outFile				string
	expectedFile		string
	defaultColumnNumber		int
}{
	{"testdata/inputFileFake.tsv", "testdata/findReplaceFake.tsv", "testdata/outputFileCreatedWholeFile.tsv", "testdata/outputFileExpectedWholeFile.tsv", -1},
}

func TestFindAndReplaceWholeFileScan(t *testing.T) {
	var err error
	for _, v := range findAndReplaceTestWholeFileScan{
		findAndReplace(v.inFile, v.findAndReplaceFile, v.outFile, v.defaultColumnNumber)
		if !fileio.AreEqual(v.outFile, v.expectedFile){
			t.Errorf("Error in findAndReplace")
		}else{
			err = os.Remove(v.outFile)
			exception.PanicOnErr(err)
		}
	}
}