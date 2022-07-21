package browser

import (
	"github.com/vertgenlab/gonomics/fileio"
	"testing"
)

var multiFaVisualizerTests = []struct {
	inFile       string
	outFile      string
	start        int
	end          int
	noMask       bool
	lineLength   int
	expectedFile string
	endOfAlignment bool
}{ //here the infile is the output region plus or minus 500 bp.
	{"testdata/chr9.112139.115209.fa", "testdata/actual.chr9.112639.114709.txt", 500, 2672, true, 200, "testdata/actual.chr9.112639.114709.txt", false},
	{"testdata/chr9.112139.115209.fa", "testdata/actual.Mask.chr9.112639.114709.txt", 500, 2672, false, 200, "testdata/maskExpected.chr9.112639.114709.txt", false},
	{"testdata/chr9.112139.115209.fa", "testdata/actual.ShortLine.chr9.112639.114709.txt", 500, 2672, false, 100, "testdata/expectedShortLine.chr9.112639.114709.txt", false},
}

func TestMultiFaVisualizer(t *testing.T) {
	for _, test := range multiFaVisualizerTests {
		MultiFaVisualizer(test.inFile, test.outFile, test.start, test.end, test.noMask, test.lineLength, test.endOfAlignment)
		if !fileio.AreEqual(test.outFile, test.expectedFile) {
			t.Errorf("Error in TestMultiFaVisualizer: expected and actual do not match.")
		}
		fileio.EasyRemove(test.outFile)
	}
}
