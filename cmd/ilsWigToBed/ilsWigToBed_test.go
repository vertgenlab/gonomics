package main

import (
	"testing"

	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/fileio"
)

var IlsWigToBedTests = []struct {
	PostProbWig    string
	ChromName      string
	ChromSizesFile string
	OutDir         string
	Expected       string
}{
	{PostProbWig: "testdata/ilsPostProbs.txt",
		ChromName:      "blah",
		ChromSizesFile: "testdata/ilsChromSizes.chrom.sizes",
		OutDir:         "testdata/ilsRecon.pfa",
		Expected:       "testdata/ilsRecon_Expected.pfa",
	},
}

func TestIlsWigToBed(t *testing.T) {
	var s IlsWigToBedSettings
	var observed []bed.Bed
	var expected []bed.Bed
	for _, v := range IlsWigToBedTests {

		s = IlsWigToBedSettings{
			PostProbWig: v.PostProbWig,
			ChromName:   v.ChromName,
			ChromSizes:  v.ChromSizesFile,
			OutPath:     v.OutDir,
		}

		IlsWigToBed(s)
		observed = bed.Read(v.OutDir)
		expected = bed.Read(v.Expected)

		if !bed.AllAreEqual(observed, expected) {
			t.Errorf("Error in ilsReconstructSeq.")
		} else {
			fileio.EasyRemove(v.OutDir)
		}
	}
}
