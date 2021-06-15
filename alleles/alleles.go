// Package alleles provides functions for counting the bases present in alignment files sam and calling variants based on those counts.
package alleles

import (
	"github.com/vertgenlab/gonomics/dna"
)

// Pileup links a genomic coordinate to an allele count for a single sample. Map structure: map[Coordinate] = AlleleCount
type Pileup map[Coordinate]*AlleleCount

// The Coordinate struct encodes the a genomic position on a linear or graph genome. In the case of a graph genome the Chr field corresponds to the node name.
type Coordinate struct {
	Chr string // or node
	Pos int
}

// Allele wraps a genomic coordinate to the allele count.
type Allele struct {
	Coordinate
	AlleleCount
}

// The AlleleCount struct holds count information for all alleles seen in an alignment file and whether those counts were on the forward or reverse read.
type AlleleCount struct {
	Ref    dna.Base
	Counts int32
	BaseAF int32
	BaseCF int32
	BaseGF int32
	BaseTF int32
	BaseAR int32
	BaseCR int32
	BaseGR int32
	BaseTR int32
	Indel  []Indel
}

// The Indel struct stores unique indels encountered in an alignment. Note that the ref and alt fields
// follow the VCF standard where the first element in the slice is the base prior to the indel
type Indel struct {
	Ref    []dna.Base
	Alt    []dna.Base
	CountF int32
	CountR int32
}