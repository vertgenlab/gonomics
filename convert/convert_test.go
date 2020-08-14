package convert

import(
	"testing"
	"github.com/vertgenlab/gonomics/vcf"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/dna"
)


var seqA []dna.Base = dna.StringToBases("--TTTC--ATGAATAA")
var seqB []dna.Base = dna.StringToBases("CCATTCCAA--CAGA-")
var inputFa []*fasta.Fasta = []*fasta.Fasta{{"eggplant", seqA}, {"squash", seqB}}
var var1 vcf.Vcf = vcf.Vcf{Chr: "chr1", Pos: 1, Id: ".", Ref: "T", Alt: "A", Qual: 100.0, Filter: "PASS", Info: ".", Format: ".", Notes: "."}
var var2 vcf.Vcf = vcf.Vcf{Chr: "chr1", Pos: 4, Id: ".", Ref: "C", Alt: "CCA", Qual: 100.0, Filter: "PASS", Info: ".", Format: ".", Notes: "."}
var var3 vcf.Vcf = vcf.Vcf{Chr: "chr1", Pos: 5, Id: ".", Ref: "ATG", Alt: "A", Qual: 100.0, Filter: "PASS", Info: ".", Format: ".", Notes: "."}
var var4 vcf.Vcf = vcf.Vcf{Chr: "chr1", Pos: 8, Id: ".", Ref: "A", Alt: "C", Qual: 100.0, Filter: "PASS", Info: ".", Format: ".", Notes: "."}
var var5 vcf.Vcf = vcf.Vcf{Chr: "chr1", Pos: 10, Id: ".", Ref: "T", Alt: "G", Qual: 100.0, Filter: "PASS", Info: ".", Format: ".", Notes: "."}
var expected []*vcf.Vcf = []*vcf.Vcf{&var1, &var2, &var3, &var4, &var5}

func TestPairwiseFaToVcf(t *testing.T) {
	input := PairwiseFaToVcf(inputFa, "chr1")
	if !vcf.AllEqual(input, expected) {
		t.Errorf("Pairwise VCF results do not match.")
	}
}