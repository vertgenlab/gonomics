package main

import (
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"os"
	"testing"
)

var CompareBedDistanceOnNameTests = []struct {
	inputBed   string
	genomeBed  string
	outputBed	string
	expectedBed   string
}{
	{"testdata/inputBed1.bed", "testdata/genomeBed.bed", "testdata/outputBed1.bed", "testdata/expectedBed1.bed"},
}

func TestCompareBedDistanceOnName(t *testing.T) {
	var err error
	for _, v := range CompareBedDistanceOnNameTests {
		compareBedOnName(v.inputBed, v.genomeBed, v.outputBed)
		if !fileio.AreEqual(v.expectedBed, v.outputBed) {
			t.Errorf("Error in compareBedDistanceBasedOnName")
		}else {
			err = os.Remove(v.outputBed)
			exception.PanicOnErr(err)
		}
	}
}



