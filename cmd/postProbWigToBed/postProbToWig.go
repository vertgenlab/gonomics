package main

import (
	"testing"

	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/fileio"
)

var PostProbToWigTests = []struct {
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

func TestPostProbToWig(t *testing.T) {
	var s PostProbToWigSettings
	var observed []bed.Bed
	var expected []bed.Bed
	for _, v := range PostProbToWigTests {

		s = PostProbToWigSettings{
			PostProbWig: v.PostProbWig,
			ChromName:   v.ChromName,
			ChromSizes:  v.ChromSizesFile,
			OutPath:     v.OutDir,
		}

		PostProbToWig(s)
		observed = bed.Read(v.OutDir)
		expected = bed.Read(v.Expected)

		if !bed.AllAreEqual(observed, expected) {
			t.Errorf("Error in ilsReconstructSeq.")
		} else {
			fileio.EasyRemove(v.OutDir)
		}
	}
}
