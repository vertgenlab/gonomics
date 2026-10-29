package reconstruct

import (
	"fmt"
	"os"
	"testing"

	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/fileio"
)

var postProbIntergenicToFaTests = []struct {
	intergenicRegions string
	inFa              string
	outPrefix         string
	expectedDir       string
	inChrom           string
}{
	{"testdata/postProbIntergenicToFa_input_1.bed", "testdata/postProbIntergenicToFa_ref_1.fasta", "testdata/postProbIntergenicToFa_test_out_1", "testdata/postProbIntergenicToFa_test_expected_1.fasta", "chr1"},
}

func TestPostProbIntergenicToFa(t *testing.T) {
	var err error
	var list_out [][]fasta.Fasta
	var out []fasta.Fasta
	var inRegions []bed.Bed
	var inFa []fasta.Fasta
	var outPath string

	for _, v := range postProbIntergenicToFaTests {
		inRegions = bed.Read(v.intergenicRegions)
		inFa = fasta.Read(v.inFa)
		list_out, out = PostProbExtractIntergenic(inRegions, inFa, v.inChrom)
		for idx, rec := range list_out {
			outPath = fmt.Sprintf("%s_%d.fasta", v.outPrefix, idx)
			fasta.Write(outPath, rec)
		}
		outPath = fmt.Sprintf("%s.fasta", v.outPrefix)
		fasta.Write(outPath, out)
		if !fileio.AreEqual(outPath, v.expectedDir) {
			t.Errorf("Error in maf MultifaToMaf")
		} else {
			err = os.Remove(outPath)
			exception.PanicOnErr(err)
			for idx := range list_out {
				outPath = fmt.Sprintf("%s_%d.fasta", v.outPrefix, idx)
				err = os.Remove(outPath)
				exception.PanicOnErr(err)
			}
		}
		exception.PanicOnErr(err)
	}
}
