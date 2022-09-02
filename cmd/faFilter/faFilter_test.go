package main

import (
	"github.com/vertgenlab/gonomics/common"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/fasta"
	"os"
	"testing"
)

var FaFilterTests = []struct {
	inputFile    string
	outputFile   string
	expectedFile string
	binExpected  string
	name         string
	notName      string
	refPositions bool
	start        int
	end          int
	minSize      int
	binFasta     int
	keepWhole    bool
}{
	{"testdata/minSizeTest.fa", "testdata/minSizeOutput.fa", "testdata/minSizeExpected.fa", "testdata/binExpected.fa", "", "", false, 0, -1, 10, 2, true},
}
var binOutputs = []string{
	"testdata/minSizeOutput.0.fa",
	"testdata/minSizeOutput.1.fa",
}

func TestFaFilter(t *testing.T) {
	var err error
	for _, v := range FaFilterTests {
		faFilter(v.inputFile, v.outputFile, v.name, v.notName, v.refPositions, v.start, v.end, v.minSize, v.binFasta, v.keepWhole)
		records := fasta.Read(v.outputFile)
		expected := fasta.Read(v.expectedFile)
		if !fasta.AllAreEqual(records, expected) {
			t.Errorf("Error in faFilter.")
		}
		err = os.Remove(v.outputFile)
		if err != nil {
			common.ExitIfError(err)
		}
		b := fasta.Read(v.binExpected)
		for r := range binOutputs {
			a := fasta.Read(binOutputs[r])
			if r == 0 {
				if fasta.IsEqual(a[0], b[0]) {
					os.Remove(binOutputs[r])
				} else {
					t.Error("1: Binning Error in faFilter.")
				}
			} else {
				if dna.CompareSeqsCaseSensitive(a[0].Seq, b[1].Seq) == 0 {
					os.Remove(binOutputs[r])
				} else {
					t.Error("2: Binning Error in faFilter.")
				}
			}
		}
	}
}
