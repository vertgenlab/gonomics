package vcf

import (
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/gtf"
)

type Variant struct {
	Vcf
	RefId   string // e.g. NC_000023.10, LRG_199, NG_012232.1, NM_004006.2, LRG-199t1, NR_002196.1, NP_003997.1, etc.
	Gene    string
	CDNAPos int
	AAPos   int
	AARef   []dna.AminoAcid
	AAAlt   []dna.AminoAcid
}

func VcfToVariant(v *Vcf, gene *gtf.Gene) *Variant {
	answer := new(Variant)

	return answer
}

func VariantToHGVS(variant *Variant) string {

}