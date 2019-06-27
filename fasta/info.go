package fasta

import (
	"github.com/vertgenlab/gonomics/dna"
	"log"
)

func Extract(f *Fasta, start int64, end int64, name string) *Fasta {
	var ans Fasta
	if start < 0 || end < 0 || start > int64(len(f.Seq)) || end > int64(len(f.Seq)) || start > end {
		log.Fatalf("Invalid start and end given to fasta.Extract. Start: %d, End: %d, ChromLength: %d.", start, end, len(f.Seq))
	}
	ans.Seq = dna.Extract(f.Seq, start, end)
	ans.Name = name
	return &ans
}
