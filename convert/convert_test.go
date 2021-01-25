package convert

import (
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/vcf"
	"strings"
	"testing"
)

var seqA, _ = dna.StringToBases("--TTTC--ATGAATAA")
var seqB, _ = dna.StringToBases("CCATTCCAA--CAGA-")
var inputFa []*fasta.Fasta = []*fasta.Fasta{{"eggplant", seqA}, {"squash", seqB}}
var var1 vcf.Vcf = vcf.Vcf{Chr: "chr1", Pos: 1, Id: ".", Ref: "T", Alt: strings.Split("A", ","), Qual: 100.0, Filter: "PASS", Info: "."}
var var2 vcf.Vcf = vcf.Vcf{Chr: "chr1", Pos: 4, Id: ".", Ref: "C", Alt: strings.Split("CCA", ","), Qual: 100.0, Filter: "PASS", Info: "."}
var var3 vcf.Vcf = vcf.Vcf{Chr: "chr1", Pos: 5, Id: ".", Ref: "ATG", Alt: strings.Split("A", ","), Qual: 100.0, Filter: "PASS", Info: "."}
var var4 vcf.Vcf = vcf.Vcf{Chr: "chr1", Pos: 8, Id: ".", Ref: "A", Alt: strings.Split("C", ","), Qual: 100.0, Filter: "PASS", Info: "."}
var var5 vcf.Vcf = vcf.Vcf{Chr: "chr1", Pos: 10, Id: ".", Ref: "T", Alt: strings.Split("G", ","), Qual: 100.0, Filter: "PASS", Info: "."}
var expected []*vcf.Vcf = []*vcf.Vcf{&var1, &var2, &var3, &var4, &var5}

func TestPairwiseFaToVcf(t *testing.T) {
	input := PairwiseFaToVcf(inputFa, "chr1")
	if !vcf.AllEqual(input, expected) {
		t.Errorf("Pairwise VCF results do not match.")
	}
}
