package fasta

import (
	"github.com/vertgenlab/gonomics/dna"
)

// Remove fasta record with index i from slice of fasta.
func Remove(slice []Fasta, i int) []Fasta {
	return append(slice[:i], slice[i+1:]...)
}

// RemoveGaps from all fasta records in a slice.
func RemoveGaps(records []Fasta) []Fasta {
	for i := range records {
		records[i].Seq = dna.RemoveGaps(records[i].Seq)
	}
	return records
}

// ReverseComplement the sequence in a fasta record.
func ReverseComplement(record Fasta) {
	dna.ReverseComplement(record.Seq)
}

// ReverseComplementAll sequences in a slice of fasta records.
func ReverseComplementAll(records []Fasta) {
	for i := range records {
		ReverseComplement(records[i])
	}
}

// ToUpper converts all bases in a fasta sequence to uppercase.
func ToUpper(fa Fasta) {
	dna.AllToUpper(fa.Seq)
}

// AllToUpper converts all bases to uppercase in all sequences
// in a slice of fasta records.
func AllToUpper(records []Fasta) {
	for i := range records {
		ToUpper(records[i])
	}
}
