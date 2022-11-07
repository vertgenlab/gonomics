package main

import (
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"os"
	"testing"
)

var BedGraphToWigTests = []struct {
	InFile       string
	ChromFile    string
	OutFile      string
	Missing      float64
	ExpectedFile string
}{
	{"testdata/test.bedGraph", "testdata/ref.chrom.sizes", "testdata/bedGraphToWig.tmp.wig", -10, "testdata/bedGraphToWig.expected.wig"},
}

func TestBedGraphToWig(t *testing.T) {
	var err error
	for _, v := range BedGraphToWigTests {
		bedGraphToWig(v.InFile, v.ChromFile, v.OutFile, v.Missing)
		if !fileio.AreEqual(v.OutFile, v.ExpectedFile) {
			t.Errorf("Error in beGraphToWig. Output is not as expected.")
		} else {
			err = os.Remove(v.OutFile)
			exception.PanicOnErr(err)
		}
	}
}
