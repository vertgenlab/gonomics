package convert

import (
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/vcf"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/exception"
	"os"
	"strings"
	"testing"
)

var seqA []dna.Base = dna.StringToBases("--TTTC--ATGAATAATCA")
var seqB []dna.Base = dna.StringToBases("CCATTCCAA--CAGAATNA")
var inputFa []*fasta.Fasta = []*fasta.Fasta{{"eggplant", seqA}, {"squash", seqB}}
var var1 vcf.Vcf = vcf.Vcf{Chr: "chr1", Pos: 1, Id: ".", Ref: "T", Alt: strings.Split("A", ","), Qual: 100.0, Filter: "PASS", Info: ".", Format: []string{"."}}
var var2 vcf.Vcf = vcf.Vcf{Chr: "chr1", Pos: 4, Id: ".", Ref: "C", Alt: strings.Split("CCA", ","), Qual: 100.0, Filter: "PASS", Info: ".", Format: []string{"."}}
var var3 vcf.Vcf = vcf.Vcf{Chr: "chr1", Pos: 5, Id: ".", Ref: "ATG", Alt: strings.Split("A", ","), Qual: 100.0, Filter: "PASS", Info: ".", Format: []string{"."}}
var var4 vcf.Vcf = vcf.Vcf{Chr: "chr1", Pos: 8, Id: ".", Ref: "A", Alt: strings.Split("C", ","), Qual: 100.0, Filter: "PASS", Info: ".", Format: []string{"."}}
var var5 vcf.Vcf = vcf.Vcf{Chr: "chr1", Pos: 10, Id: ".", Ref: "T", Alt: strings.Split("G", ","), Qual: 100.0, Filter: "PASS", Info: ".", Format: []string{"."}}
var var6 vcf.Vcf = vcf.Vcf{Chr: "chr1", Pos: 14, Id: ".", Ref: "C", Alt: strings.Split("N", ","), Qual: 100.0, Filter: "PASS", Info: ".", Format: []string{"."}}
var expected []*vcf.Vcf = []*vcf.Vcf{&var1, &var2, &var3, &var4, &var5}
var expectedSubOnly []*vcf.Vcf = []*vcf.Vcf{&var1, &var4, &var5}
var expectedRetainN []*vcf.Vcf = []*vcf.Vcf{&var1, &var2, &var3, &var4, &var5, &var6}

func TestPairwiseFaToVcf(t *testing.T) {//this test is for the default settings.
	var err error
	out := fileio.EasyCreate("tmp.txt")
	PairwiseFaToVcf(inputFa, "chr1", *out, false, false)
	input := vcf.Read("tmp.txt")
	if !vcf.AllEqual(input, expected) {
		t.Errorf("Pairwise VCF results do not match.")
	}
	err = os.Remove("tmp.txt")
	if err != nil {
		exception.PanicOnErr(err)
	}
}

func TestPairwiseFaToVcfRetainN(t *testing.T) {
	var err error
	out := fileio.EasyCreate("tmpRetainN.txt")
	PairwiseFaToVcf(inputFa, "chr1", *out, false, true)
	input := vcf.Read("tmpRetainN.txt")
	if !vcf.AllEqual(input, expectedRetainN) {
		t.Errorf("Pairwise VCF results do not match in retainN test.")
	}
	err = os.Remove("tmpRetainN.txt")
	if err != nil {
		exception.PanicOnErr(err)
	}
}

func TestPairwiseFaToVcfSubstitutionsOnly(t *testing.T) {
	var err error
	out := fileio.EasyCreate("tmpSub.txt")
	PairwiseFaToVcf(inputFa, "chr1", *out, true, false)
	input := vcf.Read("tmpSub.txt")
	if !vcf.AllEqual(input, expectedSubOnly) {
		t.Errorf("Pairwise VCF results do not match in subsitutionsOnly test.")
	}
	err = os.Remove("tmpSub.txt")
	if err != nil {
		exception.PanicOnErr(err)
	}
}

