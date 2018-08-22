package qDna

import (
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/fasta"
)

func FromBase(b dna.Base) *QBase {
	var curr QBase
	switch b {
	case dna.A:
		curr = QBase{A: 1, C: 0, G: 0, T: 0}
	case dna.C:
		curr = QBase{A: 0, C: 1, G: 0, T: 0}
	case dna.G:
		curr = QBase{A: 0, C: 0, G: 1, T: 0}
	case dna.T:
		curr = QBase{A: 0, C: 0, G: 0, T: 1}
	case dna.N, dna.Gap:
		curr = QBase{A: 0.25, C: 0.25, G: 0.25, T: 0.25}
	default:
		curr = QBase{A: 0.25, C: 0.25, G: 0.25, T: 0.25}
	}
	return &curr
}

func FromDna(in []dna.Base) []*QBase {
	answer := make([]*QBase, len(in))
	for i, _ := range in {
		answer[i] = FromBase(in[i])
	}
	return answer
}

//Convert Bases straight to QFrag --added by eric
func FromDnaToQFrag(in []dna.Base, s string) *QFrag {
	loc := Location{Assembly: "", Chr: s, Start: 0, End: 0}
	answer := QFrag{Seq: FromDna(in), From: []*Location{&loc}, Fwd: nil, Rev: nil}
	return &answer

}

func FromFasta(in *fasta.Fasta) *QFrag {
	loc := Location{Assembly: "", Chr: in.Name, Start: 0, End: 0}
	answer := QFrag{Seq: FromDna(in.Seq), From: []*Location{&loc}, Fwd: nil, Rev: nil}
	return &answer
}

func FromFastaSlice(in []*fasta.Fasta) []*QFrag {
	answer := make([]*QFrag, len(in))
	for i, _ := range in {
		answer[i] = FromFasta(in[i])
	}
	return answer
}

func mostLikelyBase(b *QBase) dna.Base {
	if b.A >= b.C && b.A >= b.G && b.A >= b.T {
		return dna.A
	} else if b.C >= b.G && b.C >= b.T {
		return dna.C
	} else if b.G >= b.T {
		return dna.G
	} else {
		return dna.T
	}
}

func mostLikelySeq(in []*QBase) []dna.Base {
	answer := make([]dna.Base, len(in))
	for i, _ := range in {
		answer[i] = mostLikelyBase(in[i])
	}
	return answer
}

func ToFasta(in *QFrag) *fasta.Fasta {
	var name string = ""
	if in.From != nil && in.From[0] != nil {
		name = in.From[0].Chr
	}
	return &fasta.Fasta{Name: name, Seq: mostLikelySeq(in.Seq)}
}



// TODO: This should be improved, but is only used
// for testing right now.  Does not follow links
func toFastaList(in []*QFrag) []*fasta.Fasta {
	answer := make([]*fasta.Fasta, len(in))
	for i, _ := range in {
		answer[i] = ToFasta(in[i])
	}
	return answer
}
