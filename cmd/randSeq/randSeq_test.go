package main

import (
	"github.com/vertgenlab/gonomics/common"
	"github.com/vertgenlab/gonomics/fasta"
	"os"
	"testing"
)

var RandSeqTests = []struct {
	outFile      string
	expectedFile string
	GC           float64
	numSeq       int
	lenSeq       int
	randSeed     bool
	setSeed      int64
}{
	{"testdata/test.fa", "testdata/expected.fa", 0.41, 10, 500, false, 10},
	{"testdata/test.fa", "testdata/expectedHighGC.fa", 0.60, 10, 500, false, 10},
	{"testdata/test.fa", "testdata/expectedShort.fa", 0.41, 10, 20, false, 10},
	{"testdata/test.fa", "testdata/expectedNumSeq.fa", 0.41, 3, 500, false, 10},
}

func TestRandSeq(t *testing.T) {
	var err error
	for _, v := range RandSeqTests {
		randSeq(v.outFile, v.GC, v.numSeq, v.lenSeq, v.randSeed, v.setSeed)
		records := fasta.Read(v.outFile)
		expected := fasta.Read(v.expectedFile)
		if !fasta.AllAreEqual(records, expected) {
			t.Errorf("Error in randSeq.")
		}
		err = os.Remove(v.outFile)
		if err != nil {
			common.ExitIfError(err)
		}
	}
}
