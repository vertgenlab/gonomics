package bed

import (
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/fasta"
)

func UngappedRegionsFromFa(fa *fasta.Fasta) []*Bed {
	var answer []*Bed
	var inRegion bool = false
	var startIndex = 0
	for index, _ := range fa.Seq {
		if dna.DefineBase(fa.Seq[index]) && inRegion == false {
			inRegion = true
			startIndex = index
		} else if !(dna.DefineBase(fa.Seq[index])) && inRegion == true {
			answer = append(answer, &Bed{Chrom: fa.Name, ChromStart: int64(startIndex), ChromEnd: int64(index)})
			inRegion = false
		}

	}
	if inRegion == true {
		answer = append(answer, &Bed{Chrom: fa.Name, ChromStart: int64(startIndex), ChromEnd: int64(len(fa.Seq))})
	}
	return answer
}

func UngappedRegionsAllFromFa(records []*fasta.Fasta) []*Bed {
	var answer []*Bed
	for idx, _ := range records {
		answer = append(answer, UngappedRegionsFromFa(records[idx])...)
	}
	return answer
}

func TotalSize(b []*Bed) int64 {
	var ans, curLen int64
	for i := 0; i < len(b); i++ {
		curLen = b[i].ChromEnd - b[i].ChromStart
		ans += curLen
	}
	return ans
}