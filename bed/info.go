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

//Uses bed regions to fragment fasta and returns fasta sequences contained in regions
//Since beds import fasta here, this function is renting until it finds a more permanant home
//Could also go in CMD folder
func EditFastaByRegion(fa *fasta.Fasta, beds []*Bed) *fasta.Fasta {
	var ans *fasta.Fasta = &fasta.Fasta{Name: fa.Name, Seq: make([]dna.Base, 0, len(fa.Seq))}
	for i, b := range beds {
		ans.Seq = append(ans.Seq, fa.Seq[b.ChromStart:b.ChromEnd]...)
		//adds 100n in between bed regions
		if i < len(beds)-2 {
			ans.Seq = append(ans.Seq, dna.CreateAllNs(100)...)
		}
	}
	return ans
}
