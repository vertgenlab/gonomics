package qDna

import (
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/fastq"
	"log"
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
	for i := range in {
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

//Convert Bases straight to QFrag --added by eric
func QFragCoord(in []dna.Base, s string, start int64, end int64) *QFrag {
	loc := Location{Assembly: "", Chr: s, Start: start, End: end}
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
	for i := range in {
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
	for i := range in {
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
	for i := range in {
		answer[i] = ToFasta(in[i])
	}
	return answer
}

func FromFastq(fq *fastq.Fastq) []*QBase {
	answer := FromBaseCalls(fq.Seq, fastq.ErrorRate(fq.Qual))
	return answer
}

func FromBaseCall(b dna.Base, err float64) *QBase {
	var curr QBase
	var e float64
	var probA float64
	var probC float64
	var probG float64
	var probT float64

	switch b {

	case dna.A:
		probA = 1 - err
		e = err / 3
		curr = QBase{A: probA, C: e, G: e, T: e}
	case dna.C:
		probC = 1 - err
		e = err / 3
		curr = QBase{A: e, C: probC, G: e, T: e}
	case dna.G:
		probG = 1 - err
		e = err / 3
		curr = QBase{A: e, C: e, G: probG, T: e}
	case dna.T:
		probT = 1 - err
		e = err / 3
		curr = QBase{A: e, C: e, G: e, T: probT}
	case dna.N:
		curr = QBase{A: 0.25, C: 0.25, G: 0.25, T: 0.25}
	default:
		log.Fatalf("Error, fastq records should not contain a gap")
	}
	return &curr
}

func FromBaseCalls(in []dna.Base, err []float32) []*QBase {
	//answer := make([]*QBase, len(in))
	if len(in) != len(err) {
		log.Fatalf("Number of bases do not match the number of quality scores")
	}
	var answer []*QBase = make([]*QBase, len(in))
	for i := 0; i < len(in); i++ {
		answer[i] = FromBaseCall(in[i], float64(err[i]))
	}
	return answer
}
