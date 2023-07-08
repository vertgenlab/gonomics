package vcf

import (
	"strings"

	"github.com/vertgenlab/gonomics/dna"
)

// FixAllVcf runs FixVcf on each element in a slice of vcf structs.
// Along with the slice of vcfs, it also needs the reference genome
// in the form of chromosome names mapping to DNA sequences.
func FixAllVcf(query []Vcf, ref map[string][]dna.Base) {
	for i := 0; i < len(query); i++ {
		FixVcf(query[i], ref)
	}
}

// FixVcf "fixes" vcf records that have a dash for a deletion,
// which does not conform to the current VCF file specs, but is often
// seen in the output of different programs.
// The function takes the vcf record to be fixed and the reference
// genome as a map of chromosome name to DNA sequence.
func FixVcf(query Vcf, ref map[string][]dna.Base) Vcf {
	return fixDash(query, ref)
}

// According to VCF specs, the alt field should never contain a "-"
// to represent a deletion. This function will reformat the vcf record
// according to VCF specs.
func fixDash(query Vcf, ref map[string][]dna.Base) Vcf {
	for i := 0; i < len(query.Alt); i++ {
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
