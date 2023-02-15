package bed

import (
	"fmt"
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

func mergeKeepLowNameAndScore(records []Bed) []Bed {
	var dist int
	var outList []Bed = make([]Bed, 0)
	SortByCoord(records)
	var curr = records[0]
	for i := 1; i < len(records); i++ {
		if Overlap(curr, records[i]) {
			if records[i].Score < curr.Score {
				curr = records[i]
			}
		} else if curr.Chrom == records[i].Chrom {
			dist = records[i].ChromStart - curr.ChromEnd
			if (curr.Score + dist) < records[i].Score { //if the distance to the record on the left plus its hidden score
				// is less than the hidden score stores for the right record, the right record will be reassigned the
				//hidden that is equal to the dist to the left record plus it's hidden value and the gene name stored on the left record
				records[i].Score = curr.Score + dist
				records[i].Name = curr.Name
			} else if (records[i].Score + dist) < curr.Score {
				curr.Score = records[i].Score + dist
				curr.Name = records[i].Name
			}
			outList = append(outList, curr)
			curr = records[i]
		} else {
			outList = append(outList, curr)
			curr = records[i]
		}
	}
	return outList
}

func removeRecordsOnEmptyChrom(records []Bed, genome map[string]chromInfo.ChromInfo) []Bed {
	var out = make([]Bed, 0)
	//var exist bool
	var v chromInfo.ChromInfo

	for i := range records {
		v, _ = genome[records[i].Chrom]
		if v.Size != 0 {
			//fmt.Printf("Check Passed: %s", v.Name)
			out = append(out, records[i])
		}
		//if exist {
		//	out = append(out, records[i])
		//}
	}
	return out
}

// FillSpaceNoHiddenValue accepts a bed and reference genome as a map[string]chromInfo.ChromInfo and returns a bed
// that assigns each genomic position to the nearest feature in the input bed using the Name field in the input bed.
// Absolute start positions of the original input bed entries are stored in the Score field of the output bed entries.
func FillSpaceNoHiddenValue(records []Bed, genome map[string]chromInfo.ChromInfo) []Bed {
	records = mergeKeepLowNameAndScore(records)
	records = removeRecordsOnEmptyChrom(records, genome)
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

//helper fucntion of FillSpaceHiddenValue. Define the 3dMidpoint between bed entries A and B (with hidden values
//Ha and Hb) as (A.ChromEnd-Ha + B.ChromStart+Hb) / 2 s.t. A.ChromEnd < B.ChromStart.
// if 3dMidpoint(A, B) < A.ChromEnd, then entry A can never have an entry in the output bed of FillSpaceHiddenValue,
// and should be removed.
func removeBedsWithNoTerritory(records []Bed) ([]Bed, bool) {
	var answer []Bed
	var threeDmidpoint int
	var violation = false

	fmt.Printf("LenRecords: %v.\n", len(records))
	for i := 1; i < len(records); i++ {
		if records[i-1].Chrom == records[i].Chrom {
			threeDmidpoint = (records[i-1].ChromEnd - records[i-1].Score + records[i].ChromStart + records[i].Score) / 2
			if threeDmidpoint < records[i-1].ChromEnd {
				violation = true
			} else {
				answer = append(answer, records[i-1])
			}
		} else {
			answer = append(answer, records[i-1])
		}
	}
	answer = append(answer, records[len(records)-1])

	return answer, violation
}

// FillSpaceHiddenValue accepts a slice of Bed structs and a reference genome as a map[string]chromInfo.ChromInfo and returns a
// slice of Bed structs that assigns each genomic position to the nearest feature in 3D space from the input bed, using
// the input bed scores to represent "hidden values", or the distance from that position to its nearest TSS in 3D space.
func FillSpaceHiddenValue(records []Bed, genome map[string]chromInfo.ChromInfo) []Bed {
	records = mergeKeepLowNameAndScore(records)
	records = removeRecordsOnEmptyChrom(records, genome)
	var violation bool = true
	for violation {
		records, violation = removeBedsWithNoTerritory(records)
	}
	var answer = make([]Bed, 0)
	var threeDMidpoint int
	var currAnswer Bed = Bed{Chrom: records[0].Chrom, ChromStart: 0, ChromEnd: records[0].ChromEnd, Name: records[0].Name, FieldsInitialized: 4}

	for i := 1; i < len(records); i++ {
		if records[i].Chrom != currAnswer.Chrom {
			currAnswer.ChromEnd = genome[records[i-1].Chrom].Size
			if currAnswer.ChromEnd-currAnswer.ChromStart < 0 {
				fmt.Printf("Here's the record in question.\n%v\n", currAnswer)
				fmt.Printf("Records[i-2]: %v.\n", records[i-2])
				fmt.Printf("Records[i-1]: %v.\n", records[i-1])
				fmt.Printf("Records[i]: %v.\n", records[i])
				fmt.Printf("Midpoint: %v.\n", threeDMidpoint)
				log.Fatalf("Died on new chrom.")
			}
			answer = append(answer, currAnswer)
			currAnswer.Chrom = records[i].Chrom
			currAnswer.ChromStart = 0
		} else if currAnswer.Name == records[i].Name {
			currAnswer.ChromStart = numbers.Min(currAnswer.ChromStart, records[i].ChromStart)
			currAnswer.ChromEnd = numbers.Max(currAnswer.ChromEnd, records[i].ChromEnd)
			currAnswer.Score = numbers.Min(currAnswer.Score, records[i].Score)
		} else {
			fmt.Printf("CurrAnswer: %v.\n", currAnswer)
			threeDMidpoint = (currAnswer.ChromEnd - records[i-1].Score + records[i].ChromStart + records[i].Score) / 2
			currAnswer.ChromEnd = threeDMidpoint + 1
			currAnswer.Name = records[i-1].Name
			if currAnswer.ChromEnd-currAnswer.ChromStart < 0 {
				fmt.Printf("Here's the record in question.\n%v\n", currAnswer)
				fmt.Printf("Records[i-2]: %v.\n", records[i-2])
				fmt.Printf("Records[i-1]: %v.\n", records[i-1])
				fmt.Printf("Records[i]: %v.\n", records[i])
				fmt.Printf("Midpoint: %v.\n", threeDMidpoint)
				log.Fatalf("Died in loop.")
			}
			answer = append(answer, currAnswer)
			currAnswer.ChromStart = threeDMidpoint + 1
		}
		currAnswer.ChromEnd = records[i].ChromEnd
		currAnswer.Name = records[i].Name
	}
	currAnswer.ChromEnd = genome[records[len(records)-1].Chrom].Size
	if currAnswer.ChromEnd-currAnswer.ChromStart < 0 {
		fmt.Printf("Here's the record in question.\n%v\n", currAnswer)
		log.Fatalf("Died after loop.")
	}
	answer = append(answer, currAnswer)
	var mergedAnswer = make([]Bed, 0)
	currAnswer = answer[0]
	for i := 1; i < len(answer); i++ {
		if answer[i-1].Name == answer[i].Name {
			currAnswer.ChromEnd = answer[i].ChromEnd
		} else {
			mergedAnswer = append(mergedAnswer, currAnswer)
			currAnswer = answer[i]
		}
	}
	mergedAnswer = append(mergedAnswer, currAnswer)
	return mergedAnswer
}
