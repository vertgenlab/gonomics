package fastq

import (
	"github.com/vertgenlab/gonomics/dna"
)

func ReverseComplement(record *Fastq) {
	dna.ReverseComplement(record.Seq)
	reverseQualRecord(record.Qual)
}

func ReverseComplementAll(records []*Fastq) {
	for idx, _ := range records {
		ReverseComplement(records[idx])
	}
}

func reverseQualRecord(qualScore []rune) {
	for i, j := 0, len(qualScore)-1; i <= j; i, j = i+1, j-1 {
		qualScore[i], qualScore[j] = qualScore[j], qualScore[i]
	}
}