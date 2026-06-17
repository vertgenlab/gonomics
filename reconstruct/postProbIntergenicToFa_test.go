package reconstruct

import (
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/exception"
 	"os"
	"testing"
)


var postProbIntergenicToFaTests = []struct {
	intergenicRegions string
	inFa	string
	outDir	string
	expectedDir	string
	inChrom string
}{
	{"testdata/postProbIntergenicToFa_input_1.bed", "testdata/postProbIntergenicToFa_ref_1.fasta", "testdata/postProbIntergenicToFa_test_out_1.fasta", "testdata/postProbIntergenicToFa_test_expected_1.fasta", "chr1",},
}

func TestPostProbIntergenicToFa(t *testing.T) {
	var err error
	var out []fasta.Fasta
	var inRegions []bed.Bed
	var inFa []fasta.Fasta

	for _, v := range postProbIntergenicToFaTests {
		inRegions = bed.Read(v.intergenicRegions)
		inFa = fasta.Read(v.inFa)
		out = PostProb_intergenic_regions_extract(inRegions, inFa, v.inChrom)
		fasta.Write(v.outDir, out)
		if !fileio.AreEqual(v.outDir, v.expectedDir) {
			t.Errorf("Error in maf MultifaToMaf")
		} else {
			err = os.Remove(v.outDir)
		}
		exception.PanicOnErr(err)
	}
}