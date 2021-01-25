package fasta

import (
	"github.com/vertgenlab/gonomics/dna"
	"log"
	"strings"
)

func Extract(f *Fasta, start int, end int, name string) *Fasta {
	var ans Fasta
	if start < 0 || end < 0 || start > len(f.Seq) || end > len(f.Seq) || start > end {
		log.Fatalf("Invalid start and end given to fasta.Extract. Start: %d, End: %d, ChromLength: %d.", start, end, len(f.Seq))
	}
	ans.Seq = dna.Extract(f.Seq, start, end)
	ans.Name = name
	return &ans
}

func CountBase(fa *Fasta, b dna.Base) int64 {
	var answer int64 = 0
	for i := 0; i < len(fa.Seq); i++ {
		if fa.Seq[i] == b {
			answer++
		}
	}
	return answer
}

func FindFaIndex(subFa []*Fasta, n string) int {
	for i := 0; i < len(subFa); i++ {
		if subFa[i].Name == n {
			return i
		}
	}
	return -1
}
func IsFasta(filename string) bool {
	if strings.HasSuffix(filename, ".fasta") || strings.HasSuffix(filename, ".fa") || strings.HasSuffix(filename, ".fasta.gz") || strings.HasSuffix(filename, ".fa.gz") {
		return true
	} else {
		return false
	}
}
