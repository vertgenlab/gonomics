package vcf

import (
	"fmt"
	"io"
	"strings"

	"github.com/vertgenlab/gonomics/dna"
)

// String implements the fmt.Stringer interface for easy printing of Vcf with the fmt package.
func (v Vcf) String() string {
	return fmt.Sprintf("%s\t%v\t%s\t%s\t%s\t%v\t%s\t%s\t%s\t%s\n", v.Chr, v.Pos, v.Id, v.Ref, strings.Join(v.Alt, ","), v.Qual, v.Filter, v.Info, v.Format, SamplesToString(v.Samples))
}

// String implements the fmt.Stringer interface for easy printing of Sample with the fmt package.
func (s Sample) String() string {
	return sampleToString(s)
}

// GetChrom returns the chromosome that the variant is on.
func (v Vcf) GetChrom() string {
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
// and adjust accordingly.

// GetChromStart returns the start position of a variant. Since VCF is 1-base, GetChromStart of a SNP variant will return v.Pos - 1.
// VCF format defines the first base of an indel as the base prior to the change, so for indels, GetChromStart will simply return v.Pos
func (v Vcf) GetChromStart() int {
	refBases := dna.StringToBases(v.Ref)
	if len(refBases) == 1 {
		return v.Pos - 1
	} else {
		return v.Pos
	}
}

// GetChromEnd returns the end position of a variant. Since VCF is 1-base, GetChromEnd of a SNP variant will rerun v.Pos.
// For indels, GetChromEnd will return v.Pos + the length of the indel - 1.
func (v Vcf) GetChromEnd() int {
	refBases := dna.StringToBases(v.Ref)
	if len(refBases) == 1 {
		return v.Pos
	} else {
		return v.Pos + len(refBases) - 1
	}
}

// UpdateCoord modifies the position data in a VCF struct based on input values.
func (v Vcf) UpdateCoord(c string, start int, end int) interface{} {
	v.Chr = c
	v.Pos = start + 1
	return v
}

// WriteToFileHandle writes a VCF struct to an io.Writer.
func (v Vcf) WriteToFileHandle(file io.Writer) {
	WriteVcf(file, v)
}
