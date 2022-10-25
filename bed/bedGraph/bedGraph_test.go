package bedGraph

import (
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"os"
	"testing"
)

var b1 BedGraph = BedGraph{Chrom: "chr1", ChromStart: 100, ChromEnd: 200, DataValue: 40.5}
var b2 BedGraph = BedGraph{Chrom: "chr2", ChromStart: 400, ChromEnd: 900, DataValue: 10.2}
var b3 BedGraph = BedGraph{Chrom: "chr3", ChromStart: 945, ChromEnd: 1000, DataValue: -3.4}
var bedGraphs []BedGraph = []BedGraph{b1, b2, b3}

var readWriteTests = []struct {
	testFileName string
	expectedFilename string
	data     []BedGraph
}{
	{"testdata/tmp.bedGraph", "testdata/bedGraphFileTest.bedGraph", bedGraphs},
}

func TestWrite(t *testing.T) {
	for _, test := range readWriteTests {
		Write(test.testFileName, test.data)
		if !fileio.AreEqual(test.testFileName, test.expectedFilename) {
			t.Errorf("The %s file was not written correctly.", test.testFileName)
		} else {
			err := os.Remove(test.testFileName)
			exception.PanicOnErr(err)
		}
	}
}

func TestRead(t *testing.T) {
	for _, test := range readWriteTests {
		read := Read(test.expectedFilename)
		for i := range read {
			if read[i].Chrom != bedGraphs[i].Chrom {
				t.Errorf("Error in reading bedGraphs. Chrom did not match.")
			}
			if read[i].ChromStart != bedGraphs[i].ChromStart {
				t.Errorf("Error in reading bedGraphs. ChromStart did not match.")
			}
			if read[i].ChromEnd != bedGraphs[i].ChromEnd {
				t.Errorf("Error in reading bedGraphs. ChromEnd did not match.")
			}
			if read[i].DataValue != bedGraphs[i].DataValue {
				t.Errorf("Error in reading bedGraphs. DataValue did not match.")
			}
		}
	}
}
