package vcf

import (
	"github.com/vertgenlab/gonomics/dna"
	"strings"
)

func FixAllVcf(query []Vcf, ref map[string][]dna.Base) {
	for i := 0; i < len(query); i++ {
		FixVcf(query[i], ref)
	}
}

func FixVcf(query Vcf, ref map[string][]dna.Base) Vcf {
	return fixDash(query, ref)
}

// According to VCF specs, the alt field should never contain a "-"
// to represent a deletion. This function will reformat the vcf record
// according to VCF specs
func fixDash(query Vcf, ref map[string][]dna.Base) Vcf {
	for i := 0; i < len(query.Alt); i++ { //TODO: (Riley) right now we're just checking if any of the alts have dashes. Not sure if this is the best way to handle polyallelic sites for fix.
		if strings.Compare(query.Alt[i], "-") == 0 {
			// query.Pos is -2, one for zero base, and one for previous base
			prevBase := dna.BaseToString(ref[query.Chr][query.Pos-2])
			query.Pos--
			query.Ref = prevBase + query.Ref
			query.Alt[i] = prevBase
		}
		return query
	}
	if strings.Compare(query.Ref, "-") == 0 {
		// query.Pos is -2, one for zero base, and one for previous base
		prevBase := dna.BaseToString(ref[query.Chr][query.Pos-2])
		query.Pos--
		query.Ref = prevBase
		for i := 0; i < len(query.Alt); i++ {
			query.Alt[i] = prevBase + query.Alt[i]
		}
	}
	return query
}
