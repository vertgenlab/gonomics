package main

import (
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

	//TODO: check line by line of out.txt, then easy remove

	//fileio.EasyRemove(out)
}
