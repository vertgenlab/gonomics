package fasta

import (
	"fmt"
	"github.com/vertgenlab/gonomics/dna"
	"log"
	"strings"
)

func AppendToName(record *Fasta, addition string) {
	record.Name = fmt.Sprintf("%s%s", record.Name, addition)
}

func AppendToNameAll(records []*Fasta, addition string) {
	for idx, _ := range records {
		AppendToName(records[idx], addition)
	}
}

func Remove(slice []*Fasta, i int) []*Fasta {
	if i < 0 || i >= len(slice) {
		log.Fatalf("Index out of range")
	}
	return append(slice[:i], slice[i+1:]...)
}

func RemoveGaps(records []*Fasta) []*Fasta {
	for i := 0 ; i < len(records); i++ {
		records[i].Seq = dna.RemoveGaps(records[i].Seq)
	}
	return records
}

func FilterName(records []*Fasta, name string) []*Fasta {
	for i := 0; i < len(records); {
		fmt.Printf("i: %d. len: %d\n", i, len(records))
		if strings.Compare(records[i].Name, name) != 0 {
			records = Remove(records, i)
		} else {
			i++
		}
	}
	return records
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

func ToUpper(fa *Fasta) {
	dna.AllToUpper(fa.Seq)
}

func AllToUpper(records []*Fasta) {
	for i := 0; i < len(records); i++ {
		ToUpper(records[i])
	}
}
