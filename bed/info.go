package bed

import (
	"fmt"

	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/fasta"
)

// UngappedRegionsFromFa: finds all regions outside gaps in a given fasta record
func UngappedRegionsFromFa(fa fasta.Fasta) []Bed {
	var answer []Bed
	var inRegion bool = false
	var startIndex, index int = 0, 0
	for index = range fa.Seq {
		if dna.DefineBase(fa.Seq[index]) && inRegion == false {
			inRegion = true
			startIndex = index
		} else if !(dna.DefineBase(fa.Seq[index])) && inRegion == true {
			answer = append(answer, Bed{Chrom: fa.Name, ChromStart: startIndex, ChromEnd: index, Name: fmt.Sprintf("%s_%d_%d", fa.Name, startIndex, index), FieldsInitialized: 4})
			inRegion = false
		}

	}
	if inRegion == true {
		answer = append(answer, Bed{Chrom: fa.Name, ChromStart: startIndex, ChromEnd: len(fa.Seq), Name: fmt.Sprintf("%s_%d_%d", fa.Name, startIndex, len(fa.Seq)), FieldsInitialized: 4})
	}
	return answer
}

// UngappedRegionsAllFromFa: Finds ungapped regions or bases that do not contain Ns. Returns a slice of bed records.
func UngappedRegionsAllFromFa(records []fasta.Fasta) []Bed {
	var answer []Bed
	var idx int = 0
	for idx = range records {
		answer = append(answer, UngappedRegionsFromFa(records[idx])...)
	}
	return answer
}

// TotalSize gives back to total region covered by bed entry.
func TotalSize(b []Bed) int {
	var ans, curLen int
	for i := 0; i < len(b); i++ {
		curLen = b[i].ChromEnd - b[i].ChromStart
		ans += curLen
	}
	return ans
}

// IsNonOverlapping returns true if any elements in a Bed slice overlap another element in the same slice. False otherwise.
// b must be presorted with SortByCoord. Verbose > 0 reveals debug prints.
func IsSelfOverlapping(b []Bed, verbose int) bool {
	for i := 0; i < len(b)-1; i++ {
		if Overlap(b[i], b[i+1]) {
			if verbose > 0 {
				fmt.Printf("first bed: %v\n", b[i])
				fmt.Printf("second bed: %v\n", b[i+1])
			}
			return true
		}
	}
	return false
}
