package vcf

import "github.com/vertgenlab/gonomics/dna"

// Current methods satisfy requirements for the following interfaces:
// bed.BedLike

func (v *Vcf) GetChrom() string {
	return v.Chr
}

// to conform with bed standards the startpos will be zero base and the endpos will be 1 base
// Example:
// bed reg  |----|
// 1-base 1  2  3  4
// 0-base 0  1  2  3
// seq    A  C  G  T
// If we have a vcf of 1 ACG -> A
// The return would be
// START = 1
// END = 3

// vcf is read in as 1-base so subtract 1 from v.pos
// for indels, vcf records the startpos as the base prior to the change
// to find the region actually being changed we need to check if it is indel
// and adjust accordingly
func (v *Vcf) GetChromStart() int {
	refBases := dna.StringToBases(v.Ref)
	if len(refBases) == 1 {
		return int(v.Pos - 1)
	} else {
		return int(v.Pos)
	}
}

func (v *Vcf) GetChromEnd() int {
	refBases := dna.StringToBases(v.Ref)
	if len(refBases) == 1 {
		return int(v.Pos)
	} else {
		return int(v.Pos) + len(refBases) - 1
	}
}
