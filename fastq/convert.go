package fastq

import (
	"github.com/vertgenlab/gonomics/dna"
)

// ReverseComplement produces the reverse complement sequence and corresponding Quals for an input Fastq record.
func ReverseComplement(record Fastq) {
	dna.ReverseComplement(record.Seq)
	ReverseQualUint8Record(record.Qual)
}

// ReverseComplementAll reverses each sequence and quality score for an input slice of Fastq records.
func ReverseComplementAll(records []Fastq) {
	for idx := range records {
		ReverseComplement(records[idx])
	}
}

func reverseQualRecord(qualScore []rune) {
	for i, j := 0, len(qualScore)-1; i <= j; i, j = i+1, j-1 {
		qualScore[i], qualScore[j] = qualScore[j], qualScore[i]
	}
}
