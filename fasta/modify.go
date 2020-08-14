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
	for i := 0; i < len(records); i++ {
		records[i].Seq = dna.RemoveGaps(records[i].Seq)
	}
	return records
}

func RefPosToAlnPos(record *Fasta, RefPos int) int {
	var AlnPos int = 0
	var refCounter int = 0

	for t := 0; refCounter < RefPos; t++ {
		AlnPos++
		if t == len(record.Seq) {
			log.Fatalf("Ran out of chromosome.")
		} else if record.Seq[t] != dna.Gap {
			refCounter++
		}
	}
	return AlnPos
}

//AlnPosToRefPos returns the reference position associated with a given AlnPos for an input Fasta. If the AlnPos corresponds to a gap, it gives the preceeding reference position.
//0 based.
func AlnPosToRefPos(record *Fasta, AlnPos int) int {
	var RefPos int = 0
	for t := 0; t < AlnPos; t++ {
		if t == len(record.Seq) {
			log.Fatalf("Ran out of chromosome.")
		} else if record.Seq[t] != dna.Gap {
			RefPos++
		}
	}
	return RefPos
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

func CopySubset(records []*Fasta, start int, end int) []*Fasta {
	c := make([]*Fasta, len(records))
	length := end - start
	for i := 0; i < len(records); i++ {
		c[i] = &Fasta{Name: records[i].Name}
		c[i].Seq = make([]dna.Base, length)
		copy(c[i].Seq, records[i].Seq[start:end])
	}
	return c
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

func GetChromIndex(records []*Fasta, name string) int {
	for i := 0; i < len(records); i++ {
		if records[i].Name == name {
			return i
		}
	}
	log.Fatalf("Chromosome name not found in fasta file.\n")
	return -1
}

//In a multiple alignment block, removes any entries comprised only of gaps.
func RemoveMissingMult(records []*Fasta) []*Fasta {
	var answer []*Fasta
	var missing bool = true

	for i := 0; i < len(records); i++ {
		missing = true
		for j := 0; j < len(records[i].Seq) && missing; j++ {
			if records[i].Seq[j] != dna.Gap {
				missing = false
			}
		}
		if !missing {
			answer = append(answer, records[i])
		}
	}
	return answer
}

//returns alignment columns with no gaps or lowercase letters
func DistColumn(records []*Fasta) []*Fasta {
	var subFa = make([]*Fasta, len(records))
	for i := 0; i < len(records); i++ {
		subFa[i] = &Fasta{Name: records[i].Name, Seq: make([]dna.Base, 0)}
	}

	for i := 0; i < len(records[0].Seq); i++ {
		currentBase := records[0].Seq[i]
		allValid := true
		if !(dna.IsLower(currentBase) || currentBase == dna.Gap) {
			for j := 1; j < len(records); j++ {
				if records[j].Seq[i] == dna.Gap || dna.IsLower(records[j].Seq[i]) {
					allValid = false
				}
			}
		} else {
			allValid = false
		}
		if allValid {
			for k := 0; k < len(records); k++ {
				subFa[k].Seq = append(subFa[k].Seq, records[k].Seq[i])
			}
		}
	}
	return subFa
}

//This function takes in a multiFa alignment block and returns only the columns that contain segregating sites.
func SegregatingSites(aln []*Fasta) []*Fasta {
	var answer = make([]*Fasta, len(aln))
	for i := 0; i < len(aln); i++ {
		answer[i] = &Fasta{Name: aln[i].Name, Seq: make([]dna.Base, 0)}
	}
	var current dna.Base
	var isSegregating bool
	for i := 0; i < len(aln[0].Seq); i++ {
		current = aln[0].Seq[i]
		isSegregating = false
		for j := 1; j < len(aln); j++ {
			if aln[j].Seq[i] != current {
				isSegregating = true
			}
		}
		if isSegregating {
			for k := 0; k < len(aln); k++ {
				answer[k].Seq = append(answer[k].Seq, aln[k].Seq[i])
			}
		}
	}
	return answer
}

func ChangePrefix(records []*Fasta, prefix string) {
	for idx := 0; idx < len(records); idx++ {
		records[idx].Name = fmt.Sprintf("%s_%d", prefix, idx)
	}
}
