package fasta

import (
	"github.com/vertgenlab/gonomics/dna"
	"strings"
)

// Copy makes a copy of the fasta struct and the associated slice of dna so that
// modifying one will not also modify the other
func Copy(a Fasta) Fasta {
	var ans Fasta
	ans = Fasta{a.Name, make([]dna.Base, len(a.Seq))}
	copy(ans.Seq, a.Seq)
	return ans
}

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

// TrimName retains the first space delimited field of a fasta name.
func TrimName(fa Fasta) Fasta {
	fields := strings.Split(fa.Name, " ")
	fa.Name = fields[0]
	return fa
}

// AllToUpper converts all bases to uppercase in all sequences
// in a slice of fasta records.
func AllToUpper(records []Fasta) {
	for i := range records {
		ToUpper(records[i])
	}
}
