package fasta

import (
	"fmt"
	"github.com/vertgenlab/gonomics/dna"
)

func AppendToName(record *Fasta, addition string) {
	record.Name = fmt.Sprintf("%s%s", record.Name, addition)
}

func AppendToNameAll(records []Fasta, addition string) {
	for idx, _ := range records {
		AppendToName(&records[idx], addition)
	}
}

func ReverseComplement(record *Fasta) {
	dna.ReverseComplement(record.Seq)
}

func ReverseComplementAll(records []Fasta) {
	for idx, _ := range records {
		ReverseComplement(&records[idx])
	}
}
