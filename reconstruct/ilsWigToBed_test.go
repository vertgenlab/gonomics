package reconstruct

import (
	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/wig"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/exception"
	"testing"
	"strconv"
	"os"
)

var ilsWigToBedTests = []struct {
	postProbWigFile string
	chromSizes	string
	outDirPrefix	string
	expectedPrefix	string
}{
	{"testdata/postProbToWig_Out_1.wig", "testdata/postProbs.sizes", "testdata/ilsWigToBed_Out_v", "testdata/ilsWigToBed_Expected_v"},
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