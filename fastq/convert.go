package fastq

import (
	"github.com/vertgenlab/gonomics/dna"
	//"github.com/vertgenlab/gonomics/qDna"
	//"log"
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


/*
func FromFastq(fq *Fastq) []*qDna.QBase {
	answer := FromDna(fq.Seq, ErrorRate(fq.Qual))
	return answer
}*/