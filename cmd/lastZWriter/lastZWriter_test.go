package main

import (
	"github.com/vertgenlab/gonomics/fileio"
	"strings"
	"testing"
)

var lastZ = "lastZInstall"
var pairwise = "../../lastZWriter/testdata"
var speciesListFile = pairwise + "/speciesList.txt"
var refListFile = pairwise + "/refList.txt"
var allDists = pairwise + "/allDistsAll.txt"
var out = "out.txt"

func TestMakeArray(t *testing.T) {
	MakeArray(lastZ, pairwise, speciesListFile, refListFile, allDists, out, true, "")
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
