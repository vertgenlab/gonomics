package bed

import (
	"github.com/vertgenlab/gonomics/numbers"
	"log"
)

// TrimSlice trims each element in a slice of bed entries by an input left and right amount.
func TrimSlice(b []Bed, trimLeft int, trimRight int) {
	for i := range b {
		b[i] = Trim(b[i], trimLeft, trimRight)
	}
}

// Trim shortens bed entries on the left and right side by an input-specified number of bases. These values must not exceed the length of the bed entry and must be non-negative.
func Trim(b Bed, trimLeft int, trimRight int) Bed {
	if trimLeft < 0 || trimRight < 0 {
		log.Fatalf("Error in bed/Trim. Must trim bed values by a value greater or equal to zero.")
	}
	b.ChromStart = b.ChromStart + trimLeft
	b.ChromEnd = b.ChromEnd - trimRight
	if b.ChromStart >= b.ChromEnd {
		log.Fatalf("Error in Trim. Attempted to remove too much from bed entry. Please select a lower trim value or exclude the bed entry as position %v\t%v.\n", b.Chrom, b.ChromStart)
	}
	return b
}

// AllToMidpoint edits all input Bed structs in an input slice so that its coordinates correspond to its midpoint.
func AllToMidpoint(b []Bed) {
	for i := range b {
		b[i] = ToMidpoint(b[i])
	}
}

// ToMidpoint edits an input Bed struct so that its coordinates correspond to its midpoint.
func ToMidpoint(b Bed) Bed {
	midpoint := (b.ChromStart + b.ChromEnd) / 2
	b.ChromStart = midpoint
	b.ChromEnd = midpoint + 1
	return b
}

// ToTss edits an input Bed struct so that its coordinates corresponds to the start position, strand-sensitive.
func ToTss(b Bed) Bed {
	switch b.Strand {
	case Positive:
		b.ChromEnd = b.ChromStart + 1
	case Negative:
		b.ChromStart = b.ChromEnd - 1
	default:
		log.Fatalf("Input bed must have an annotated positive or negative strand to trim to Tss.")
	}
	return b
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
func MergeHighMem(records []Bed, mergeAdjacent int, keepAllNames bool) []Bed {
	var outList []Bed
	var minDist int
	var err error
	if len(records) == 0 {
		return records //empty and nil slices are returned as is.
	}
	SortByCoord(records)
	var currentMax = records[0]

	for i := 1; i < len(records); i++ {
		minDist, err = MinimumDistance(currentMax, records[i])
		if Overlap(currentMax, records[i]) || (minDist <= mergeAdjacent && err == nil) {
			if records[i].Score > currentMax.Score {
				currentMax.Score = records[i].Score
			}
			currentMax.ChromEnd = numbers.Max(records[i].ChromEnd, currentMax.ChromEnd)
			if keepAllNames && records[i].Name != "" {
				if currentMax.Name != "" {
					currentMax.Name = currentMax.Name + "," + records[i].Name
				} else {
					currentMax.Name = records[i].Name
				}
			}
		} else {
			outList = append(outList, currentMax)
			currentMax = records[i]
		}
	}
	outList = append(outList, currentMax)
	return outList
}
