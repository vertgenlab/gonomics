package fasta

import (
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/bed"
)

func UngappedRegions(fa *Fasta) ([]bed.Bed) {
	var answer []bed.Bed
	var inRegion bool = false
	var startIndex = 0
 	for index, _ := range fa.Seq {
 		if dna.DefineBase(fa.Seq[index]) && inRegion == false {
 			inRegion = true
 			startIndex = index
 		} else if !(dna.DefineBase(fa.Seq[index])) && inRegion == true {
				answer = append(answer, bed.Bed{Name: fa.Name, Start:startIndex, End:index})
				inRegion = false
			}	
		 
	}
	if inRegion == true {
		answer = append(answer, bed.Bed{Name: fa.Name, Start:startIndex, End:len(fa.Seq)})
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
