package sam

import (
	"github.com/vertgenlab/gonomics/cigar"
	"github.com/vertgenlab/gonomics/dna"
)

// Equal returns true if the two input Aln structs are identical.
func Equal(a Aln, b Aln) bool {
	if a.QName != b.QName {
		return false
	}
	if a.Flag != b.Flag {
		return false
	}
	if a.RName != b.RName {
		return false
	}
	if a.Pos != b.Pos {
		return false
	}
	if a.MapQ != b.MapQ {
		return false
	}
	if cigar.ToString(a.Cigar) != cigar.ToString(b.Cigar) {
		return false
	}
	if a.RNext != b.RNext {
		return false
	}
	if a.PNext != b.PNext {
		return false
	}
	if a.TLen != b.TLen {
		return false
	}
	if dna.CompareSeqsIgnoreCase(a.Seq, b.Seq) != 0 {
		return false
	}
	if a.Qual != b.Qual {
		return false
	}
	if a.Extra != b.Extra {
		return false
	}
	return true
}
