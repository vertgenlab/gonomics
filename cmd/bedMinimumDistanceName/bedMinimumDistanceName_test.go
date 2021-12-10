package main

import (
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"os"
	"testing"
)

var BedMinimumDistanceNameTests = []struct {
	inputBed    string
	genomeBed   string
	outputBed   string
	expectedBed string
}{
	{"testdata/inputBed1.bed", "testdata/genomeBed1.bed", "testdata/outputBed1.bed", "testdata/expectedBed1.bed"},
	{"testdata/inputBed2.bed", "testdata/genomeBed2.bed", "testdata/outputBed2.bed", "testdata/expectedBed2.bed"},
	{"testdata/inputBed3.bed", "testdata/genomeBed3.bed", "testdata/outputBed3.bed", "testdata/expectedBed3.bed"},
}

func TestBedMinimumDistanceName(t *testing.T) {
	var err error
	for _, v := range BedMinimumDistanceNameTests {
		bedMinimumDistanceName(v.inputBed, v.genomeBed, v.outputBed)
		if !fileio.AreEqual(v.expectedBed, v.outputBed) {
			t.Errorf("Error in bedMinimumDistanceName")
		} else {
			err = os.Remove(v.outputBed)
			exception.PanicOnErr(err)
		}
	}
}

