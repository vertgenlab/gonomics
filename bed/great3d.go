package bed

import (
	"github.com/vertgenlab/gonomics/chromInfo"
	"github.com/vertgenlab/gonomics/numbers"
	"log"
)

// FillSpaceNoHiddenValue accepts a bed of single point tss positions and reference genome as a map[string]chromInfo.ChromInfo and returns a bed
// that assigns each genomic position to the nearest feature in the input bed using the Name field in the input bed.
// Absolute start positions of the original input bed entries are stored in the Score field of the output bed entries.
func FillSpaceNoHiddenValue(records []Bed, genome map[string]chromInfo.ChromInfo) []Bed {
	records = removeRecordsOnMissingChrom(records, genome)
	var answer = make([]Bed, 0)
	if records == nil || len(records) == 0 {
		return records // return slice as is if it is empty or nil.
	}

	var currAnswer = Bed{Chrom: records[0].Chrom, ChromStart: 0, ChromEnd: records[0].ChromEnd, Name: records[0].Name, Score: records[0].ChromStart, FieldsInitialized: 5}
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

// removeRecordsOnMissingChrom will remove any gene records that exist on a chromosome that isn't in the chrom.sizes file provided
func removeRecordsOnMissingChrom(records []Bed, genome map[string]chromInfo.ChromInfo) []Bed {
	var out = make([]Bed, 0)
	var exist bool

	for i := range records {
		_, exist = genome[records[i].Chrom]
		if exist {
			out = append(out, records[i])
		}
	}
	return out
}

// helper fucntion of FillSpaceHiddenValue. Define the 3dMidpoint between bed entries A and B (with hidden values
// Ha and Hb) as (A.ChromEnd-Ha + B.ChromStart+Hb) / 2 s.t. A.ChromEnd < B.ChromStart.
// if 3dMidpoint(A, B) < A.ChromEnd, then entry A can never have an entry in the output bed of FillSpaceHiddenValue,
// and should be removed.
func removeBedsWithNoTerritory(records []Bed) ([]Bed, bool) {
	var answer []Bed
	var threeDmidpoint int
	var violation = false
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
	records = runUntilNoNewHidden(records)
	records = removeRecordsOnMissingChrom(records, genome)
	var violation = true
	for violation {
		records, violation = removeBedsWithNoTerritory(records)
	}
	var answer = make([]Bed, 0)
	var threeDMidpoint int
	var currAnswer = Bed{Chrom: records[0].Chrom, ChromStart: 0, ChromEnd: records[0].ChromEnd, Name: records[0].Name, FieldsInitialized: 4}

	for i := 1; i < len(records); i++ {
		if records[i].Chrom != currAnswer.Chrom {
			currAnswer.ChromEnd = genome[records[i-1].Chrom].Size
			if currAnswer.ChromEnd < currAnswer.ChromStart {
				log.Print(records[i-1])
				log.Print(currAnswer)
				log.Fatalf("Died on new chrom.")
			}
			answer = append(answer, currAnswer)
			currAnswer.Chrom = records[i].Chrom
			currAnswer.ChromStart = 0
		} else if currAnswer.Name == records[i].Name && currAnswer.Chrom == records[i].Chrom {
			currAnswer.ChromStart = numbers.Min(currAnswer.ChromStart, records[i].ChromStart)
			currAnswer.ChromEnd = numbers.Max(currAnswer.ChromEnd, records[i].ChromEnd)
			currAnswer.Score = numbers.Min(currAnswer.Score, records[i].Score)
		} else {
			threeDMidpoint = (currAnswer.ChromEnd - records[i-1].Score + records[i].ChromStart + records[i].Score) / 2
			currAnswer.ChromEnd = threeDMidpoint + 1
			currAnswer.Name = records[i-1].Name
			if currAnswer.ChromEnd-currAnswer.ChromStart < 0 {
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
		log.Fatalf("Died after loop.")
	}
	answer = append(answer, currAnswer)

	return answer
}

// runUntilNoNewHidden directs our mergeKeepLowNameAndScore function to keep updating until no changes occur
func runUntilNoNewHidden(records []Bed) []Bed {
	var newHidden bool
	records, newHidden = mergeKeepLowScoreAndName(records)
	if newHidden {
		runUntilNoNewHidden(records)
	}
	return records
}

// mergeKeepLowScoreAndName is a helper function for 3d great that takes in a slice of bed structs
func mergeKeepLowScoreAndName(records []Bed) (out []Bed, newHidden bool) {
	var dist int
	newHidden = false
	var outList = make([]Bed, 0)
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
				// is less than the hidden score stored for the right record, the right record will be reassigned the
				//hidden value that is equal to the dist to the left record plus it's hidden value and the gene name stored on the left record
				newHidden = true
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
	outList = append(outList, curr)
	return outList, newHidden
}
