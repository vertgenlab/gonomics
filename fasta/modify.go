package fasta

import (
	"fmt"
	"github.com/vertgenlab/gonomics/dna"
)

func AppendToName(record *Fasta, addition string) {
	record.Name = fmt.Sprintf("%s%s", record.Name, addition)
}

func AppendToNameAll(records []*Fasta, addition string) {
	for idx, _ := range records {
		AppendToName(records[idx], addition)
	}
}

func ReverseComplement(record *Fasta) {
	dna.ReverseComplement(record.Seq)
}

func ReverseComplementAll(records []*Fasta) {
	for idx, _ := range records {
		ReverseComplement(records[idx])
	}
}

func DivideFasta(fa *Fasta, n int) []*Fasta {
	var answer []*Fasta
	leftover := len(fa.Seq) % n
	for i := 0; i < len(fa.Seq)-leftover; i += n {
		answer = append(answer, &Fasta{Name: fmt.Sprintf("%s_%d", fa.Name, i), Seq: fa.Seq[i : i+n]})
	}
	return answer
}

func DivideFastaAll(fa []*Fasta, n int) [][]*Fasta {
	var answer [][]*Fasta
	for index, _ := range fa {
		answer = append(answer, DivideFasta(fa[index], n))
	}
	return answer
}
