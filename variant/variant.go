package variant

import (
	"github.com/vertgenlab/gonomics/dna"
)

// Substitution describes a single base change.
// Pos is the 0-base position of the changed base.
// (e.g. ATG sub 1 T>C -> ACG)
type Substitution struct {
	Chr string
	Pos int
	Ref dna.Base
	Alt dna.Base
}

// Insertion describes an insertion into a sequence.
// Pos is the 0-base position of the base after the insertion.
// (e.g. ATG ins 1 C -> ACTG)
type Insertion struct {
	Chr string
	Pos int
	Seq []dna.Base
}

// Deletion describes a sequence deletion.
// Start and End are the 0-base left-closed right-open interval
// being deleted (e.g. ATG del [1,2) -> AG).
type Deletion struct {
	Chr   string
	Start int
	End   int
}

// Delins describes a combined insertion and a deletion of
// a sequence. Useful for describing complex variants.
// Start and End are the 0-base left-closed right-open interval
// being deleted (e.g. ATG del [1,3) ins C -> AC).
type Delins struct {
	Chr    string
	Start  int
	End    int
	InsSeq []dna.Base
}

// TODO break into multiple types???
// Structural is a catch-all for large variants including
// complex rearrangements, large insertions/deletions, and
// transposition of mobile elements.
type Structural struct {
}

// TODO Duplication???

// TODO Inversion???
