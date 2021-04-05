package variant

import (
	"errors"
	"github.com/vertgenlab/gonomics/dna"
)

var (
	ErrRefMatch        = errors.New("position in seq does not match expected ref")
	ErrInvalidPosition = errors.New("position of mutation was invalid")
	ErrNegPos          = errors.New("variant position plus offset is negative")
)

// Mutator is satisfied by types that implement the Mutate method, which alters
// an input []dna.Base per coordinates and sequence stored in the receiver type.
// The offsetStart/End integer inputs for the Mutate method are added to the start
// position and end position of the receiver type and may be positive or negative.
// Receiver types that only have one position (e.g. Substitution and Insertion)
// will only use the offsetStart value. The offset input is useful for cases
// where only a subset of the origin chromosome is available, operating on spliced
// sequences, or when making multiple insertions/deletions to the sequence.
//
// The Mutate method does not guarantee that the input sequence remains unchanged.
// The return of the Mutate method should be stored in the input variable.
// e.g. sequence := Insertion.Mutate(sequence, 0)
//
// Mutator is satisfied by Substitution, Insertion, Deletion, and Delins.
// Structural variation is often too complex to implement the simple Mutate
// method, therefore changing sequences per a structural variant should be done
// with a separate function. // TODO make this function and update godoc.
type Mutator interface {
	Mutate(seq []dna.Base, offsetStart int, offsetEnd int) ([]dna.Base, error)
}

// Mutate for Substitution makes the substitution directly in the input sequence.
// Note that this changes the input sequence, even if the return slice is saved to
// a different variable.
func (s Substitution) Mutate(seq []dna.Base, offsetStart int, offsetEnd int) ([]dna.Base, error) {
	offsetPos := s.Pos + offsetStart

	if offsetPos < 0 {
		return nil, ErrNegPos
	}

	if seq[offsetPos] != s.Ref {
		return nil, ErrRefMatch
	} else {
		seq[offsetPos] = s.Alt
	}
	return seq, nil
}

// Mutate for Insertion attempts to make the insertion directly in the input
// sequence. If the input []dna.Base does not have the capacity to store the
// inserted bases, a new slice is allocated to store the input sequence and insertion.
// Note that this may change the input sequence, even if the return slice is saved to
// a different variable.
func (i Insertion) Mutate(seq []dna.Base, offsetStart int, offsetEnd int) ([]dna.Base, error) {
	offsetPos := i.Pos + offsetStart

	if offsetPos < 0 {
		return nil, ErrNegPos
	}

	if offsetPos > len(seq) {
		return nil, ErrInvalidPosition
	}
	return dna.Insert(seq, offsetPos, i.Seq), nil
}

// Mutate for Deletion makes the deletion directly in the input sequence.
// Note that this changes the input sequence, even if the return slice is
// saved to a different variable.
func (d Deletion) Mutate(seq []dna.Base, offsetStart int, offsetEnd int) ([]dna.Base, error) {
	offsetStartPos := d.Start + offsetStart
	offsetEndPos := d.End + offsetEnd

	if offsetStartPos < 0 {
		return nil, ErrNegPos
	}

	if offsetEndPos > len(seq) {
		offsetEndPos = len(seq)
	}
	return dna.Delete(seq, offsetStartPos, offsetEndPos), nil
}

// Mutate for Delins attempts to make a combined deletion and insertion directly
// in the input sequence. If the input []dna.Base does not have the capacity to store
// the resulting sequence, a new slice is allocated. Note that this will change the
// input sequence, even if the return slice is saved to a different variable.
func (di Delins) Mutate(seq []dna.Base, offsetStart int, offsetEnd int) ([]dna.Base, error) {
	offsetStartPos := di.Start + offsetStart
	offsetEndPos := di.End + offsetEnd

	if offsetStartPos < 0 {
		return nil, ErrNegPos
	}

	if offsetEndPos > len(seq) {
		offsetEndPos = len(seq)
	}

	return dna.Replace(seq, offsetStartPos, offsetEndPos, di.InsSeq), nil
}
