package main

import (
	"os"
	"testing"

	"github.com/vertgenlab/gonomics/bed/bedpe"
)

var test = []struct {
	strawFile  string
	outFile    string
	chrom      string
	binSize    int
	interChrom string
	expected   string
}{
	{"testdata/in.straw", "testdata/out.bedpe", "chr1", 5000, "", "testdata/expected.out.bedpe"},
	{"testdata/in.straw", "testdata/out.interChrom.bedpe", "chr1", 5000, "chr2", "testdata/expected.out.interChrom.bedpe"},
}

func TestStrawToBedpe(t *testing.T) {
	for i := range test {
		strawToBedpe(test[i].strawFile, test[i].outFile, test[i].chrom, test[i].binSize, test[i].interChrom)
		if !bedpe.AllAreEqual(bedpe.Read(test[i].expected), bedpe.Read(test[i].outFile)) {
			t.Errorf("outFile: %s did not match expected file: %s.", test[i].outFile, test[i].expected)
		} else {
			os.Remove(test[i].outFile)
		}
	}
}
