package variant

import (
	"errors"
	"github.com/vertgenlab/gonomics/dna"
)

var (
	ErrRefMatch        error = errors.New("position in seq does not match expected ref")
	ErrInvalidPosition error = errors.New("position of mutation was invalid")
)

type Mutater interface {
	Mutate(seq []dna.Base) ([]dna.Base, error)
}

func (s Substitution) Mutate(seq []dna.Base) ([]dna.Base, error) {
	if seq[s.Pos] != s.Ref {
		return nil, ErrRefMatch
	} else {
		seq[s.Pos] = s.Alt
	}
	return seq, nil
}

func (i Insertion) Mutate(seq []dna.Base) ([]dna.Base, error) {
	if i.Pos > len(seq) {
		return nil, ErrInvalidPosition
	}
	return dna.Insert(seq, i.Pos, i.Seq), nil
}

func (d Deletion) Mutate(seq []dna.Base) ([]dna.Base, error) {
	if d.End > len(seq) {
		return nil, ErrInvalidPosition
	}
	return dna.Delete(seq, d.Start, d.End), nil
}

func (id Indel) Mutate(seq []dna.Base) ([]dna.Base, error) {
	seq, err := id.Deletion.Mutate(seq)
	if err != nil {
		return nil, err
	}
	return id.Insertion.Mutate(seq)
}
