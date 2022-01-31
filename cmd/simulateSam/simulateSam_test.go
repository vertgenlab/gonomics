package main

import (
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/sam"
	"math/rand"
	"os"
	"testing"
)

func TestSimulateSam(t *testing.T) {
	rand.Seed(1)
	simulateSam("testdata/test.fa", "testdata/actual.sam", 100, false)
	simulateSam("testdata/test.fa", "testdata/actual.bam", 100, true)

	if !fileio.AreEqual("testdata/actual.sam", "testdata/expected.sam") {
		t.Error("problem simulating sam")
	}

	bamA, _ := sam.Read("testdata/actual.bam")
	bamB, _ := sam.Read("testdata/expected.bam")

	if len(bamA) != len(bamB) {
		t.Error("problem simulating bam")
	}

	for i := range bamA {
		if !sam.Equal(bamA[i], bamB[i]) {
			t.Error("problem simulating bam")
		}
	}

	if !t.Failed() {
		err := os.Remove("testdata/actual.sam")
		exception.PanicOnErr(err)
		err = os.Remove("testdata/actual.bam")
		exception.PanicOnErr(err)
	}
}
