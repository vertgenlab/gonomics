package main

import (
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"os"
	"testing"
)

func TestLongReadLibStats(t *testing.T) {
	var err error
	longReadLibStats("testdata/in.fq", "testdata/out.stats.txt", "testdata/out.sizes.txt")
	if !fileio.AreEqual("testdata/out.stats.txt", "testdata/exp.stats.txt") {
		t.Errorf("Error in LongReadLibStats: library stats don't match expected\n")
	} else {
		err = os.Remove("testdata/out.stats.txt")
		exception.PanicOnErr(err)
	}

	if !fileio.AreEqual("testdata/out.sizes.txt", "testdata/exp.sizes.txt") {
		t.Errorf("Error in LongReadStats: read lengths don't match expected\n")
	} else {
		err = os.Remove("testdata/out.sizes.txt")
		exception.PanicOnErr(err)
	}

}
