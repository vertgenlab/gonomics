package main

import (
	"os"
	"testing"

	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
)

var findAndReplaceColumnTests = []struct {
	inFile              string
	inFileDelim         string
	findAndReplaceFile  string
	findAndReplaceDelim string
	outFile             string
	expectedFile        string
	columnNumber        int
	ignoreColumns       bool
}{
	{"testdata/inputFileFake.tsv", "\t", "testdata/findReplaceFake.tsv", "\t", "testdata/outputFileCreatedColumn0.tsv",
		"testdata/outputFileExpectedColumn0.tsv", 0, false},
	{"testdata/inputFileFake.tsv", "\t", "testdata/findReplaceFake.tsv", "\t", "testdata/outputFileCreatedColumn1.tsv",
		"testdata/outputFileExpectedColumn1.tsv", 1, false},
	{"testdata/inputFileFake.tsv", "\t", "testdata/findReplaceFake.tsv", "\t", "testdata/outputFileCreatedColumn2.tsv",
		"testdata/outputFileExpectedColumn2.tsv", 2, false},
	{"testdata/inputFileFake.tsv", "\t", "testdata/findReplaceFake.tsv", "\t", "testdata/outputFileCreatedWholeFile.tsv",
		"testdata/outputFileExpectedWholeFile.tsv", -1, false},
	{"testdata/inputOne.txt", "\t", "testdata/findReplaceOne.txt", "\t", "testdata/temp.txt",
		"testdata/expectedOne.txt", -1, true},
}

func TestFindAndReplaceColumnSpecific(t *testing.T) {
	var err error
	for _, v := range findAndReplaceColumnTests {
		findAndReplace(v.inFile, v.inFileDelim, v.findAndReplaceFile, v.findAndReplaceDelim, v.outFile, v.columnNumber, v.ignoreColumns)
		if !fileio.AreEqual(v.outFile, v.expectedFile) {
			t.Errorf("Error in findAndReplace")
		} else {
			err = os.Remove(v.outFile)
			exception.PanicOnErr(err)
		}
	}
}

/*
var findAndReplaceTestWholeFileScan = []struct {
	inFile              string
	findAndReplaceFile  string
	outFile             string
	expectedFile        string
	defaultColumnNumber int
}{
	{"testdata/inputFileFake.tsv", "testdata/findReplaceFake.tsv", "testdata/outputFileCreatedWholeFile.tsv", "testdata/outputFileExpectedWholeFile.tsv", -1},
}

func TestFindAndReplaceWholeFileScan(t *testing.T) {
	var err error
	for _, v := range findAndReplaceTestWholeFileScan {
		findAndReplace(v.inFile, v.findAndReplaceFile, v.outFile, v.defaultColumnNumber)
		if !fileio.AreEqual(v.outFile, v.expectedFile) {
			t.Errorf("Error in findAndReplace")
		} else {
			err = os.Remove(v.outFile)
			exception.PanicOnErr(err)
		}
	}
}*/
