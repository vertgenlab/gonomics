package main

import (
	"testing"

	"github.com/vertgenlab/gonomics/fasta/pFasta"
	"github.com/vertgenlab/gonomics/fileio"
)

var IlsReconstructTests = []struct {
	PostProbsFiles string
	ReconFiles     string
	ChromSizesFile string
	OutDir         string
	Precision      float32
	Expected       string
	TestPrecision  float32
}{
	{PostProbsFiles: "testdata/ilsPostProbs.txt",
		ReconFiles:     "testdata/ilsReconsInput.txt",
		ChromSizesFile: "testdata/ilsChromSizes.chrom.sizes",
		Precision:      0.001,
		OutDir:         "testdata/ilsRecon.pfa",
		Expected:       "testdata/ilsRecon_Expected.pfa",
		TestPrecision:  0.001,
	},
}

func TestSumIlsSeqs(t *testing.T) {
	var s SumIlsSeqsSettings
	var observed []pFasta.PFasta
	var expected []pFasta.PFasta
	for _, v := range IlsReconstructTests {

		s = SumIlsSeqsSettings{
			PostProbsFiles: v.PostProbsFiles,
			ReconFiles:     v.ReconFiles,
			ChromSizesFile: v.ChromSizesFile,
			OutDir:         v.OutDir,
			Precision:      v.Precision,
		}

		SumIlsSeqs(s)
		observed = pFasta.Read(v.OutDir)
		expected = pFasta.Read(v.Expected)

		if !pFasta.AllAreEqual(observed, expected, v.TestPrecision) {
			t.Errorf("Error in iSumIlsSeqs.")
		} else {
			fileio.EasyRemove(v.OutDir)
		}
	}
}
