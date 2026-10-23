package reconstruct

import (
	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/wig"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/exception"
	"os"
	"testing"
	"strconv"
)


var ilsWigToBedTests = []struct {
	postProbWigFile string
	chromSizes	string
	outDirPrefix	string
	expectedPrefix	string
}{
	// {"C:/Users/sarah/Documents/pDNA/pdna/hcg_postProb_wig_luria/hcg_postProb.wig", "C:/Users/sarah/Documents/pDNA/pdna/20way/hcg_postProb.sizes", "C:/Users/sarah/Documents/pDNA/pdna/20way/ilsWigToBedOut/20way_postProb_v", "C:/Users/sarah/Documents/pDNA/pdna/20way/ilsWigToBedOut/20way_postProb_v"},
	{"testdata/postProbToWig_Expected_1.wig", "testdata/ilsWigToBed_test_1.sizes", "testdata/wigToBed_test_output_v", "testdata/wigToBed_test_expected_v"},
}

func TestIlsWigToBed(t *testing.T) {
	var err error
	var out [][]bed.Bed
	var inWig map[string]wig.Wig
	for _, v := range ilsWigToBedTests {
		inWig = wig.Read(v.postProbWigFile, v.chromSizes, 0)
		out = IlsWigToBed(inWig)

		for idx, ilsBed := range out {
			if len(ilsBed) >= 0 {
				bed.Write(v.outDirPrefix+strconv.Itoa(idx) + ".bed", out[idx])
			}

			if !fileio.AreEqual(v.outDirPrefix+strconv.Itoa(idx)+".bed", v.expectedPrefix+strconv.Itoa(idx)+".bed") {
				t.Errorf("Error in reconstruct ilsWigToBed")
			} else {
				err = os.Remove(v.outDirPrefix+strconv.Itoa(idx)+".bed")
			}
			exception.PanicOnErr(err)
		}
		
	}
}