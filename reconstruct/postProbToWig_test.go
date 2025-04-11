package reconstruct

import (
	"github.com/vertgenlab/gonomics/maf"
	"github.com/vertgenlab/gonomics/wig"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/exception"
	"testing"
	"os"
)

var multiFaToMafTests = []struct {
	mafFile 	string
	postProbFile	string
	outDir	string
	expectedDir	string
}{
	{"testdata/test_maf.maf", "testdata/summed_post_prob.csv", "testdata/postProbToWig_Out_1.wig", "testdata/postProbToWig_Expected_1.wig"},
}

func TestPostProbToWig(t *testing.T) {
	var err error
	var out map[string]wig.Wig
	var inMaf []*maf.Maf
	for _, v := range multiFaToMafTests {
		inMaf = maf.Read(v.mafFile)
		out = PostProbToWig(v.postProbFile, inMaf)
		wig.Write(v.outDir, out)
		if !fileio.AreEqual(v.outDir, v.expectedDir) {
			t.Errorf("Error in maf MultifaToMaf")
		} else {
			err = os.Remove(v.outDir)
		}
		exception.PanicOnErr(err)
	} 
}