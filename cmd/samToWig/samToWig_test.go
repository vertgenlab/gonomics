package main

import (
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"os"
	"testing"
)

var SamToWigTests = []struct {
	inFile     string
	reference  string
	outFile    string
	expected   string
	fragLength int
	deletions  bool
}{
	{"testdata/test1.sam", "testdata/test.chrom.sizes", "testdata/out.test1.wig", "testdata/test1.wig", -1, false},
	{"testdata/test2.sam", "testdata/test.chrom.sizes", "testdata/out.test2.wig", "testdata/test2.wig", 30, false},
	{"testdata/test1.bam", "testdata/test.chrom.sizes", "testdata/out.test1.wig", "testdata/test1.wig", -1, false},
	{"testdata/test2.bam", "testdata/test.chrom.sizes", "testdata/out.test2.wig", "testdata/test2.wig", 30, false},
	{"testdata/test1.sam", "testdata/test.chrom.sizes", "testdata/out.test1.withDel.wig", "testdata/test1.withDel.wig", -1, true},
}

func TestSamToWig(t *testing.T) {
	var err error
	for _, v := range SamToWigTests {
		samToWig(v.inFile, v.reference, v.outFile, v.fragLength, v.deletions)
		if !fileio.AreEqual(v.expected, v.outFile) {
			t.Errorf("Error in samToWig: Exp: %s, actual: %s\n", v.expected, v.outFile)
		}
		err = os.Remove(v.outFile)
		exception.PanicOnErr(err)
	}
}
