package main

import (
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"os"
	"testing"
)

var digestGenomeTests = []struct {
	inFile       string
	outFile      string
	expectedFile string
	motif        string
}{
	{"testdata/in.fa", "testdata/out.MboI.bed", "testdata/expected.MboI.bed", "MboI"},
	{"testdata/in.fa", "testdata/out.CGCG.bed", "testdata/expected.CGCG.bed", "C^GCG"},
	{"testdata/in.fa", "testdata/out.AAGA.bed", "testdata/expected.AAGA.bed", "A^AGA"},
}

func TestDigestGenome(t *testing.T) {
	var err error
	for _, v := range digestGenomeTests {
		digestGenome(v.inFile, v.motif, v.outFile)
		if !fileio.AreEqual(v.expectedFile, v.outFile) {
			t.Errorf("Error: digestGenome output files %s and %s are not equal to one another...", v.outFile, v.expectedFile)
		} else {
			err = os.Remove(v.outFile)
			exception.PanicOnErr(err)
		}
	}
}
