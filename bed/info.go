package bed

import (
	"fmt"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/fasta"
	"log"
)

// MultiFaUngappedRegions takes in a MultiFa format alignment and returns a bed file of all ungapped regions
// for a user-specified seqName. The resulting sequence is in target coordinates, and the bed chrom field can
// be specified by the argument 'chromName'.
func MultiFaUngappedRegions(records []fasta.Fasta, chromName string, seqName string) []Bed {
	var answer = make([]Bed, 0)
	var seqNameIndex int = -1
	var inRegion bool = false
	var startRefPos, endRefPos int = 0, 0
	var lastRefPos, lastAlnPos = 0, 0
	var currAlnPos int = 0
	var seqNameFound bool = false

	//first we find the index of seqName
	for i := 0; i < len(records); i++ {
		if records[i].Name == seqName && !seqNameFound {
			seqNameIndex = i
			seqNameFound = true
		} else if records[i].Name == seqName && seqNameFound {
			log.Fatalf("Error: found the same record, %s, multiple times in input fasta.\n", seqName)
		}
	}

	if seqNameIndex == -1 {
		log.Fatalf("Error: seqName: %s, not found in records.\n", seqName)
	}

	for currAlnPos = 0; currAlnPos < len(records[0].Seq); currAlnPos++ {
		if dna.DefineBase(records[seqNameIndex].Seq[currAlnPos]) && !inRegion {
			inRegion = true
			startRefPos = fasta.AlnPosToRefPosCounter(records[0], currAlnPos, lastRefPos, lastAlnPos)
			lastRefPos, lastAlnPos = startRefPos, currAlnPos
		} else if !(dna.DefineBase(records[seqNameIndex].Seq[currAlnPos])) && inRegion {
			endRefPos = fasta.AlnPosToRefPosCounter(records[0], currAlnPos, lastRefPos, lastAlnPos)
			lastRefPos, lastAlnPos = endRefPos, currAlnPos
			answer = append(answer, Bed{Chrom: chromName, ChromStart: startRefPos, ChromEnd: endRefPos, FieldsInitialized: 3})
			inRegion = false
		}
	}

	if inRegion {
		endRefPos = fasta.AlnPosToRefPosCounter(records[0], currAlnPos, lastRefPos, lastAlnPos)
		answer = append(answer, Bed{Chrom: chromName, ChromStart: startRefPos, ChromEnd: endRefPos, FieldsInitialized: 3})
	}

	return answer
}

// UngappedRegionsFromFa finds all regions outside gaps in a given fasta record.
func UngappedRegionsFromFa(fa fasta.Fasta) []Bed {
	var answer []Bed
	var inRegion bool = false
	var startIndex, index int = 0, 0
	for index = range fa.Seq {
		if dna.DefineBase(fa.Seq[index]) && !inRegion {
			inRegion = true
			startIndex = index
		} else if !(dna.DefineBase(fa.Seq[index])) && inRegion {
			answer = append(answer, Bed{Chrom: fa.Name, ChromStart: startIndex, ChromEnd: index, Name: fmt.Sprintf("%s_%d_%d", fa.Name, startIndex, index), FieldsInitialized: 4})
			inRegion = false
		}
	}
	if inRegion {
		answer = append(answer, Bed{Chrom: fa.Name, ChromStart: startIndex, ChromEnd: len(fa.Seq), Name: fmt.Sprintf("%s_%d_%d", fa.Name, startIndex, len(fa.Seq)), FieldsInitialized: 4})
	}
	return answer
}

// UngappedRegionsAllFromFa finds ungapped regions or bases that do not contain Ns. Returns a slice of bed records.
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

// IsSelfOverlapping returns true if any elements in a Bed slice overlap another element in the same slice. False otherwise.
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

// Size returns the size of a singular bed entry as an int
func Size(b Bed) int {
	return b.ChromEnd - b.ChromStart
}
