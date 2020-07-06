package gtf

import (
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/vcf"
)

type Variant struct {
	vcf.Vcf
	RefId   string // e.g. NC_000023.10, LRG_199, NG_012232.1, NM_004006.2, LRG-199t1, NR_002196.1, NP_003997.1, etc.
	Gene    string
	CDNAPos int
	AAPos   int
	AARef   []dna.AminoAcid
	AAAlt   []dna.AminoAcid
}

func VcfToVariant(v *vcf.Vcf, genes map[string]*Gene) *Variant {
	answer := new(Variant)
	MoveAllCanonicalToZero(genes)

	return answer
}

func VariantToHGVS(variant *Variant) string {

}

