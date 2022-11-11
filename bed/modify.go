package bed

import (
	"github.com/vertgenlab/gonomics/chromInfo"
	"github.com/vertgenlab/gonomics/numbers"
	"log"
)

//Trim shortens bed entries on the left and right side by an input-specified number of bases. These values must not exceed the length of the bed entry and must be non-negative.
func Trim(b []Bed, trimLeft int, trimRight int) {
	if trimLeft < 0 || trimRight < 0 {
		log.Fatalf("Error in bed/Trim. Must trim bed values by a value greater or equal to zero.")
	}
	for i := range b {
		b[i].ChromStart = b[i].ChromStart + trimLeft
		b[i].ChromEnd = b[i].ChromEnd - trimRight
		if b[i].ChromStart >= b[i].ChromEnd {
			log.Fatalf("Error in Trim. Attempted to remove too much from bed entry. Please select a lower trim value or exclude the bed entry as position %v\t%v.\n", b[i].Chrom, b[i].ChromStart)
		}
	}
}

/*
//input must be sorted. incomplete function, does not work.
func MergeLowMem(b <- chan Bed, mergeAdjacent bool) <- chan Bed {
	var firstTime bool = true
	var out chan Bed
	var currentMax Bed

	for i := range b {
		if firstTime {
			firstTime = false
			currentMax = i
		} else {
			if Overlap(currentMax, i) || mergeAdjacent && Adjacent(currentMax, i) {
				if i.Score > currentMax.Score {
					currentMax.Score = i.Score
				}
				currentMax.ChromEnd = numbers.Max(i.ChromEnd, currentMax.ChromEnd)
			} else {
				out <- currentMax
				currentMax = i
			}
		}
	}
	return out
}
*/

// MergeHighMem retains input Bed entries that are non-overlapping with other input bed entries and merges together overlapping bed entries.
// Merged bed entries will retain the maximum score in the output.
func MergeHighMem(records []Bed, mergeAdjacent bool) []Bed {
	var outList []Bed
	if records == nil || len(records) == 0 {
		return records //empty and nil slices are returned as is.
	}
	SortByCoord(records)
	var currentMax Bed = records[0]

	for i := 1; i < len(records); i++ {
		if Overlap(currentMax, records[i]) || mergeAdjacent && Adjacent(currentMax, records[i]) {
			if records[i].Score > currentMax.Score {
				currentMax.Score = records[i].Score
			}
			currentMax.ChromEnd = numbers.Max(records[i].ChromEnd, currentMax.ChromEnd)
		} else {
			outList = append(outList, currentMax)
			currentMax = records[i]
		}
	}
	outList = append(outList, currentMax)
	return outList
}

// FillSpaceNoHiddenValue accepts a bed and reference genome as a map[string]chromInfo.ChromInfo and returns a bed
// that assigns each genomic position to the nearest feature in the input bed using the Name field in the input bed.
// Absolute start positions of the original input bed entries are stored in the Score field of the output bed entries.
func FillSpaceNoHiddenValue(records []Bed, genome map[string]chromInfo.ChromInfo) []Bed {
	SortByCoord(records)
	MergeHighMem(records, true) // possible silent errors if two features are directly adjacent. We chose to discard.
	var answer = make([]Bed, 0)

	if records == nil || len(records) == 0 {
		return records // return slice as is if it is empty or nil.
	}

	var currAnswer Bed = Bed{Chrom: records[0].Chrom, ChromStart: 0, ChromEnd: records[0].ChromEnd, Name: records[0].Name, Score: records[0].ChromStart, FieldsInitialized: 5}
	var midpoint int

	//when a position is equidistant from two features, it is assigned to the left feature.
	for i := 1; i < len(records); i++ {
		if records[i].Chrom != currAnswer.Chrom {
			currAnswer.ChromEnd = genome[currAnswer.Chrom].Size
			answer = append(answer, currAnswer)
			currAnswer.Chrom = records[i].Chrom
			currAnswer.ChromStart = 0
			currAnswer.ChromEnd = records[i].ChromEnd
			currAnswer.Name = records[i].Name
			currAnswer.Score = records[i].ChromStart
		} else {
			midpoint = (records[i].ChromStart + currAnswer.ChromEnd) / 2
			currAnswer.ChromEnd = midpoint + 1
			answer = append(answer, currAnswer)
			currAnswer.Chrom = records[i].Chrom
			currAnswer.ChromStart = midpoint + 1
			currAnswer.ChromEnd = records[i].ChromEnd
			currAnswer.Name = records[i].Name
			currAnswer.Score = records[i].ChromStart
		}
	}

	currAnswer.ChromEnd = genome[currAnswer.Chrom].Size
	answer = append(answer, currAnswer)

	return answer
}
