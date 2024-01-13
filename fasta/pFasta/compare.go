package pFasta

import (
	"github.com/vertgenlab/gonomics/dna/pDna"
	"strings"
)

// AllAreEqual returns true if two []pFasta data structures are determined
// to be equal, false otherwise.
// Two []pFasta data structures are said to be equal if: (1) they are of the same
// length. (2) Each component pFasta struct in 'a' is equal to the pFasta struct
// at the corresponding index in 'b'. See the documentation for 'Equal' for a
// complete definition of the criterion for pFasta struct equality.
// Note that pFasta equality is parameterized by a user-specified relative 'precision'.
func AllAreEqual(a []pFasta, b []pFasta, precision float32) bool {
	if len(a) != len(b) {
		return false
	}
	for i := range a {
		if !Equal(a[i], b[i], precision) {
			return false
		}
	}
	return true
}

// Equal returns true if two input pFasta structs are 'equal', false otherwise.
// Two input pFasta structs are said to be equal if they have the same name,
// have the same number of bases in their sequence, and for each base, the
// base probabilities are equal within a user-specified level of relative precision.
func Equal(a pFasta, b pFasta, precision float32) bool {
	if strings.Compare(a.Name, b.Name) != 0 {
		return false
	}
	if len(a.Seq) != len(b.Seq) {
		return false
	}
	for i := range a.Seq {
		if !pDna.EqualBase(a.Seq[i], b.Seq[i], precision) {
			return false
		}
	}
	return true
}
