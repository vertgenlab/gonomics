package main

import (
	"os"
	"testing"

	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
)

var BedDistanceFromChrEndsTests = []struct {
	inBed       string
	chromSize   string
	outBed      string
	expectedBed string
}{{"testdata/input1.bed",
	"testdata/test.chrom.sizes",
	"testdata/output1.bed",
	"testdata/expected1.bed"},
	{"testdata/input2.bed",
		"testdata/test.chrom.sizes",
		"testdata/output2.bed",
		"testdata/expected2.bed"},
}

func TestBedDistanceFromChrEnds(t *testing.T) {
	var err error
	for _, v := range BedDistanceFromChrEndsTests {
		bedDistanceFromChrEnds(v.inBed, v.chromSize, v.outBed)
		if !fileio.AreEqual(v.expectedBed, v.outBed) {
			t.Errorf("Error in bedDistanceFromChrEnds.go, expectedBed != outputBed")
		} else {
			err = os.Remove(v.outBed)
			exception.PanicOnErr(err)
		}

	}
}
