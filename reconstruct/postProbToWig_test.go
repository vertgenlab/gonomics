package reconstruct

import (
	"github.com/vertgenlab/gonomics/maf"
	"github.com/vertgenlab/gonomics/wig"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/exception"
	"testing"
	"os"
)

var postProbtoWigTests = []struct {
	mafFile 	string
	postProbFile	string
	desiredSpecies []string
	wigOutDir	string
	mapEstimOutDir	string
	expectedWigDir	string
	expectedMapDir	string
}{
	{"testdata/chr22.test.maf", "testdata/postProbToWig_post_prob_test.csv", []string{"hg38", "panPan2", "panTro6", "gorGor5"}, "testdata/postProbToWig_out.wig", "testdata/postProbToWig_map_out.wig", "testdata/postProbToWig_expected.wig", "testdata/postProbToWig_map_expected.wig"},
}

func TestPostProbToWig(t *testing.T) {
	// TODO: there is a weird print out of "diff" + list of a few numbers? probably some sort of probability 
	var err error

	var postProbOut map[string]wig.Wig
	var mapEstimOut map[string]wig.Wig
	var inMaf []*maf.Maf
	for _, v := range postProbtoWigTests {
		inMaf = maf.Read(v.mafFile)
		postProbOut, mapEstimOut  = PostProbToWig(v.postProbFile, inMaf, v.desiredSpecies)
		wig.Write(v.wigOutDir, postProbOut)
		wig.Write(v.mapEstimOutDir, mapEstimOut)
		if !fileio.AreEqual(v.wigOutDir, v.expectedWigDir) {
			t.Errorf("Error in maf MultifaToMaf wig")
		} else if !fileio.AreEqual(v.mapEstimOutDir, v.expectedMapDir) {
			t.Errorf("Error in maf MultifaToMaf map")
		} else {
			err = os.Remove(v.wigOutDir)
			err = os.Remove(v.mapEstimOutDir)
		}
		exception.PanicOnErr(err)
	} 
}