package fasta

import (
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/bed"
)

func UngappedRegions(fa *Fasta) ([]bed.Bed) {
	var answer []bed.Bed
 	for index, _ := range fa.Seq {
		if dna.DefineBase(fa.Seq[index]) {
			answer = append(answer, bed.Bed{Name: fa.Name, Start:index, End:index+1})
		}
					
	}		
	return answer
}

func UngappedRegionsAll(records []Fasta) ([]bed.Bed){
	var b []bed.Bed
        for idx, _ := range records {
               b = append(b, UngappedRegions(&records[idx])...)
	}
	return b
}
