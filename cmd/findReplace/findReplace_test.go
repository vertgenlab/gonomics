package main

import (
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"os"
	"testing"
)

var SingleLineTests = []struct {
	inputLine    string
	inputFind    []string
	inputReplace []string
	expectedLine string
}{
	{"I like to eat tacos.  I work in the CaRL Building.  Genomics is fun.", []string{"taco", "crazy", "Building.  G"}, []string{"salsa and chip", "shouldNotAppear", "building; g"}, "I like to eat salsa and chips.  I work in the CaRL building; genomics is fun."},
	{"This should not change", []string{"red", "blue", "yellow"}, []string{"yes", "no", "maybe"}, "This should not change"},
	{"elementOne", []string{"red", "elementOne", "yellow"}, []string{"yes", "elementA", "maybe"}, "elementA"},
}

func TestFindReplaceOneLine(t *testing.T) {
	var actual string
	for _, v := range SingleLineTests {
		actual = findReplaceOneLine(v.inputLine, v.inputFind, v.inputReplace)
		if actual != v.expectedLine {
			t.Errorf("Error in findReplaceOneLine.  Actual result:\"%s\" Expected result:\"%s\"\n", actual, v.expectedLine)
		}
	}
}

var FileTests = []struct {
	inputFile           string
	inputFind           string
	inputReplace        string
	outputFile          string
	expectedFile        string
	findReplaceAreFiles bool
}{
	{"testdata/inputOne.txt", "testdata/findOne.txt", "testdata/replaceOne.txt", "testdata/actual.txt", "testdata/expectedOne.txt", true},
	{"testdata/inputTwo.txt", "Blue", "Devils", "testdata/actual.txt", "testdata/expectedTwo.txt", false},
}

func TestFaFormat(t *testing.T) {
	var err error
	for _, v := range FileTests {
		findReplace(v.inputFile, v.inputFind, v.inputReplace, v.outputFile, v.findReplaceAreFiles)
		if !fileio.AreEqual(v.outputFile, v.expectedFile) {
			t.Errorf("Error in findReplace; output did not match expected.")
		} else {
			err = os.Remove(v.outputFile)
			exception.PanicOnErr(err)
		}
	}
}
