package bed

import (
	"fmt"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/fasta"
)

//UngappedRegionsFromFa: finds all regions outside gaps in a given fasta record
func UngappedRegionsFromFa(fa fasta.Fasta) []Bed {
	var answer []Bed
	var inRegion bool = false
	var startIndex, index int = 0, 0
	for index = range fa.Seq {
		if dna.DefineBase(fa.Seq[index]) && inRegion == false {
			inRegion = true
			startIndex = index
		} else if !(dna.DefineBase(fa.Seq[index])) && inRegion == true {
			answer = append(answer, Bed{Chrom: fa.Name, ChromStart: startIndex, ChromEnd: index, FieldsInitialized: 3})
			inRegion = false
		}

	}
	if inRegion == true {
		answer = append(answer, Bed{Chrom: fa.Name, ChromStart: startIndex, ChromEnd: len(fa.Seq), FieldsInitialized: 3})
	}
	return answer
}

//UngappedRegionsAllFromFa: Finds ungapped regions or bases that do not contain Ns. Returns a slice of bed records.
func UngappedRegionsAllFromFa(records []fasta.Fasta) []Bed {
	var answer []Bed
	var idx int = 0
	for idx = range records {
		answer = append(answer, UngappedRegionsFromFa(records[idx])...)
	}
	return answer
}

//TotalSize gives back to total region covered by bed entry.
func TotalSize(b []Bed) int {
	var ans, curLen int
	for i := 0; i < len(b); i++ {
		curLen = b[i].ChromEnd - b[i].ChromStart
		ans += curLen
	}
	return ans
}

//IsNonOverlapping returns true if any elements in a Bed slice overlap another element in the same slice. False otherwise.
//b must be presorted with SortByCoord.
func IsSelfOverlapping(b []Bed) bool {
	for i := 0; i < len(b)-1; i++ {
		if Overlap(b[i], b[i+1]) {
			fmt.Printf("first bed: %v\n", b[i])
			fmt.Printf("second bed: %v\n", b[i+1])
			return true
		}
	}
	return false
}

//Splits fasta regions by using bed regions and concatenate fasta sequences by filling 100 Ns in between
func MakeContigFromBed(fa *fasta.Fasta, beds []Bed) *fasta.Fasta {
	var ans *fasta.Fasta = &fasta.Fasta{Name: fa.Name, Seq: make([]dna.Base, 0)}
	for i, b := range beds {
		ans.Seq = append(ans.Seq, fa.Seq[b.ChromStart:b.ChromEnd]...)
		//adds 100n in between bed regions
		if i < len(beds)-2 {
			ans.Seq = append(ans.Seq, dna.CreateAllNs(100)...)
		}
	}
	return ans
}
