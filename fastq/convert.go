package fastq

import (
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/qDna"
	"log"
)

func ReverseComplementFastq(record *Fastq) *Fastq {
	var answer *Fastq = record
	dna.ReverseComplement(answer.Seq)
	swapQualRecord(answer.Qual)
	return answer

}

func ReverseComplementAll(records []*Fastq) []*Fastq {
	var answer []*Fastq
	for idx, _ := range records {
		answer = append(answer, ReverseComplementFastq(records[idx]))
	}
	return answer
}

func swapQual(alpha rune, beta rune) (rune, rune) {
	return beta, alpha
}

func swapQualRecord(qualScore []rune) {
	for i, j := 0, len(qualScore)-1; i <= j; i, j = i+1, j-1 {
		qualScore[i], qualScore[j] = swapQual(qualScore[i], qualScore[j])
	}
}

func FromBase(b dna.Base, err float64) *qDna.QBase {
	var curr qDna.QBase
	switch b {
	case dna.A:
		probA := 1 - err
		e := err / 3
		curr = qDna.QBase{A: probA, C: e, G: e, T: e}
	case dna.C:
		probC := 1 - err
		e := err / 3
		curr = qDna.QBase{A: e, C: probC, G: e, T: e}
	case dna.G:
		probG := 1 - err
		e := err / 3
		curr = qDna.QBase{A: e, C: e, G: probG, T: e}
	case dna.T:
		probT := 1 - err
		e := err / 3
		curr = qDna.QBase{A: e, C: e, G: e, T: probT}
	case dna.N, dna.Gap:
		curr = qDna.QBase{A: 0.25, C: 0.25, G: 0.25, T: 0.25}
	default:
		curr = qDna.QBase{A: 0.25, C: 0.25, G: 0.25, T: 0.25}
	}
	return &curr
}

func FromDna(in []dna.Base, err []float64) []*qDna.QBase {
	//answer := make([]*QBase, len(in))
	if len(in) != len(err) {
		log.Fatalf("Number of bases do not match the number of quality scores")
	}
	var answer []*qDna.QBase
	for i := 0; i < len(in); i++ {
		answer = append(answer, FromBase(in[i], err[i]))
	}
	return answer
}

func FromFastq(fq *Fastq) []*qDna.QBase {
	answer := FromDna(fq.Seq, ErrorRate(fq.Qual))
	return answer
}