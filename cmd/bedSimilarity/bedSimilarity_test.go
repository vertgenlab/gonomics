package main

import (
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"os"
	"testing"
)

var HaqerBedCompTests = []struct {
	bedA             string
	bedB             string
	outFile          string
	expFile          string
	list             string
	matrixAverage    string
	matrixComponents string
	expectedMatrix   string
}{
	{"testdata/smallAJ.bed", "testdata/largeAJ.bed", "testdata/out.twoComps.txt", "testdata/expected.twoComps.txt", "", "", "", ""},
	{"", "", "testdata/out.list.txt", "testdata/expected.list.txt", "testdata/list.txt", "testdata/out.matrixAvg.txt", "", "testdata/expected.matrixAvg.txt"},
	{"", "", "testdata/out.list.txt", "testdata/expected.list.txt", "testdata/list.txt", "", "testdata/out.matrixComp.txt", "testdata/expected.matrixComp.txt"},
}

func TestHaqerBedComps(t *testing.T) {
	var err error
	var s settings
	for _, i := range HaqerBedCompTests {
		s = settings{
			bedA:             i.bedA,
			bedB:             i.bedB,
			outFile:          i.outFile,
			list:             i.list,
			matrixComponents: i.matrixComponents,
			matrixAverage:    i.matrixAverage,
		}
		bedSimilarity(s)
		if !fileio.AreEqual(i.outFile, i.expFile) {
			t.Errorf("Error: test files %s and %s are not equal to one another...", i.outFile, i.expFile)
		} else {
			err = os.Remove(i.outFile)
			exception.PanicOnErr(err)
		}
		if i.matrixComponents != "" {
			if !fileio.AreEqual(i.matrixComponents, i.expectedMatrix) {
				t.Errorf("Error: test files %s and %s are not equal to one another...", i.matrixComponents, i.expectedMatrix)
			} else {
				err = os.Remove(i.matrixComponents)
				exception.PanicOnErr(err)
			}
		}
		if i.matrixAverage != "" {
			if !fileio.AreEqual(i.matrixAverage, i.expectedMatrix) {
				t.Errorf("Error: test files %s and %s are not equal to one another...", i.matrixAverage, i.expectedMatrix)
			} else {
				err = os.Remove(i.matrixAverage)
				exception.PanicOnErr(err)
			}
		}
	}
}
