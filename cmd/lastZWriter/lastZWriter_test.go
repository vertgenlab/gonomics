package main

import (
	"strings"
	"testing"

	"github.com/vertgenlab/gonomics/fileio"
)

var lastZ = "lastZInstall"
var pairwise = "../../lastZWriter/testdata"
var speciesListFile = pairwise + "/speciesList.txt"
var refListFile = pairwise + "/refList.txt"
var allDists = pairwise + "/allDistsAll.txt"
var out = "testdata/out.txt"

func TestMakeArray(t *testing.T) {
	MakeArray(lastZ, pairwise, speciesListFile, refListFile, allDists, out, true, "", "")
	outRecords := fileio.EasyOpen(out)
	expected := fileio.Read("testdata/expected.txt")
	lineNum := 0
	for l, done := fileio.EasyNextRealLine(outRecords); !done; l, done = fileio.EasyNextRealLine(outRecords) {
		lineNum++
		inputs := strings.Split(l, " ")
		exp := strings.Split(expected[lineNum-1], " ")
		for i := range inputs {
			match := strings.Compare(inputs[i], exp[i])
			if match != 0 {
				t.Fatalf("Output line %d: %s, did not match expected value: %s", lineNum, inputs[i], expected[i])
			}
		}
	}

	fileio.EasyRemove(out)
}

var speciesListFileSimple = pairwise + "/speciesList_simple.txt"
var refListFileSimple = pairwise + "/refList_simple.txt"
var parameters = "M=50 K=2200"
var outSimple = "testdata/out_simple.txt"
var targetModifier = "[unmask]"

func TestMakeArraySimple(t *testing.T) {
	MakeArraySimple(lastZ, pairwise, speciesListFileSimple, refListFileSimple, parameters, outSimple, targetModifier)
	outRecords := fileio.EasyOpen(outSimple)
	expected := fileio.Read("testdata/expected_simple.txt")
	lineNum := 0
	for l, done := fileio.EasyNextRealLine(outRecords); !done; l, done = fileio.EasyNextRealLine(outRecords) {
		lineNum++
		inputs := strings.Split(l, " ")
		exp := strings.Split(expected[lineNum-1], " ")
		for i := range inputs {
			match := strings.Compare(inputs[i], exp[i])
			if match != 0 {
				t.Fatalf("Error: Output line %d: %s, did not match expected value: %s", lineNum, inputs[i], expected[i])
			}
		}
	}

	fileio.EasyRemove(outSimple)
}
