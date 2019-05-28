package fasta

import (
	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/dna"
)

func UngappedRegions(fa *Fasta) []*bed.Bed {
	var answer []*bed.Bed
	var inRegion bool = false
	var startIndex = 0
	for index, _ := range fa.Seq {
		if dna.DefineBase(fa.Seq[index]) && inRegion == false {
			inRegion = true
			startIndex = index
		} else if !(dna.DefineBase(fa.Seq[index])) && inRegion == true {
			answer = append(answer, &bed.Bed{Chrom: fa.Name, ChromStart: int64(startIndex), ChromEnd: int64(index)})
			inRegion = false
		}

	}
	if inRegion == true {
		answer = append(answer, &bed.Bed{Chrom: fa.Name, ChromStart: int64(startIndex), ChromEnd: int64(len(fa.Seq))})
	}
	return answer
}

func UngappedRegionsAll(records []Fasta) []*bed.Bed {
	var answer []*bed.Bed
	for idx, _ := range records {
		answer = append(answer, UngappedRegions(&records[idx])...)
	}
	return answer
}
