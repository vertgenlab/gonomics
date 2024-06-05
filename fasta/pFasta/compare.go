package pFasta

import (
	"fmt"
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
func AllAreEqual(a []PFasta, b []PFasta, precision float32) bool {
	if len(a) != len(b) {
		return false
	}
	for i := range a {
		if !IsEqual(a[i], b[i], precision) {
			return false
		}
	}
	return true
}

// IsEqual returns true if two input pFasta structs are 'equal', false otherwise.
// Two input pFasta structs are said to be equal if they have the same name,
// have the same number of bases in their sequence, and for each base, the
// base probabilities are equal within a user-specified level of relative precision.
func IsEqual(a PFasta, b PFasta, precision float32) bool {
	if strings.Compare(a.Name, b.Name) != 0 {
		return false
	}
	if len(a.Seq) != len(b.Seq) {
		return false
	}
	for i := range a.Seq {
		if !pDna.EqualBase(a.Seq[i], b.Seq[i], precision) {
			fmt.Printf("%e\t%e\n", a.Seq[i].A, b.Seq[i].A)
			fmt.Printf("%e\t%e\n", a.Seq[i].C, b.Seq[i].C)
			fmt.Printf("%e\t%e\n", a.Seq[i].G, b.Seq[i].G)
			fmt.Printf("%e\t%e\n", a.Seq[i].T, b.Seq[i].T)
			return false
		}
	}
	return true
}

// IsValid returns false if any base position in the sequence is not valid
func IsValid(pFaIn PFasta, precision float32) bool {
	for i := range pFaIn.Seq {
		if !pDna.IsValid(pFaIn.Seq[i], precision) {
			return false
		}
	}
	return true
}
