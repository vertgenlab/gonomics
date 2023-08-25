package ontology

import (
	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/bed/bedpe"
	"github.com/vertgenlab/gonomics/chromInfo"
	"github.com/vertgenlab/gonomics/interval"
	"github.com/vertgenlab/gonomics/numbers"
	"log"
)

// Fill3dSpace takes in a set of Bedpe format contact points, single-base pair []bed.Bed representing
// transcription start sites, and a reference genome as a map[string]chromInfo.ChromInfo. This program
// returns a []bed.Bed representing the closest TSS in 3D space for each position in the reference genome,
// using the Bedpe contacts to determine the closest TSS in 3D space.
// Note that the input TSS must have 0 in score fields.
func Fill3dSpace(contacts []bedpe.BedPe, tss []bed.Bed, sizes map[string]chromInfo.ChromInfo) []bed.Bed {
	var chromInFile bool
	geneChroms := makeChromSlice(tss) //makes a list of the chroms in the gene file
	answer := make([]bed.Bed, len(tss))
	copy(answer, tss)
	var currNearestBed, currAnswerA, currAnswerB bed.Bed
	var closest1dGenesIntervals, currNearest []interval.Interval
	closest1dGene := FillSpaceNoHiddenValue(tss, sizes)
	if contacts == nil {
		return closest1dGene
	}

	for i := range closest1dGene {
		closest1dGenesIntervals = append(closest1dGenesIntervals, closest1dGene[i])
	}
	closest1dGeneTree := interval.BuildTree(closest1dGenesIntervals)

	bedpe.ContactsToMidpoints(contacts)

	// this for loop finds the nearest gene and hidden value for each bedpe foot midpoint
	for j := range contacts {
		//this is just a check to make sure that for any new chromosomes we encounter that they exist in the gene file, so we don't throw an error.
		if j == 0 {
			chromInFile = checkGeneFileForChrom(contacts[j], geneChroms)
		} else if contacts[j].A.Chrom != contacts[j-1].A.Chrom { // only check when we encounter a new chrom in the bedpe file
			chromInFile = checkGeneFileForChrom(contacts[j], geneChroms)
		}
		if !chromInFile {
			continue
		}
		currNearest = interval.Query(closest1dGeneTree, contacts[j].A, "any")
		if len(currNearest) > 1 || len(currNearest) == 0 {
			log.Fatalf("Space Filled bed should return one nearest bed entry, returned %v.", currNearest)
		}
		currNearestBed = currNearest[0].(bed.Bed)
		currAnswerA = bed.Bed{Chrom: contacts[j].A.Chrom,
			ChromStart:        contacts[j].A.ChromStart,
			ChromEnd:          contacts[j].A.ChromEnd,
			Name:              currNearestBed.Name,
			Score:             0,
			FieldsInitialized: 5}
		//we've maintained the original absolute chromStart of every gene record in the FillSpaceNoHiddenValue run so
		//that when we are calculating the distance to that start we have an easy way to find that position.
		//the actual start position contained for the record after running FillSpaceNoHiddenValue would be the midpoint between that gene and its left neighbor
		if currNearestBed.Score < contacts[j].A.ChromStart { //if the absolute start position of a gene's tss (currNearestBed.Score) is left of the bedpe midpoint
			currAnswerA.Score = contacts[j].A.ChromStart - currNearestBed.Score
		} else {
			currAnswerA.Score = currNearestBed.Score - contacts[j].A.ChromStart
		}

		currNearest = interval.Query(closest1dGeneTree, contacts[j].B, "any")
		if len(currNearest) > 1 || len(currNearest) == 0 {
			log.Fatal("Space Filled bed should return one nearest bed entry.")
		}
		currNearestBed = currNearest[0].(bed.Bed)
		currAnswerB = bed.Bed{Chrom: contacts[j].B.Chrom,
			ChromStart:        contacts[j].B.ChromStart,
			ChromEnd:          contacts[j].B.ChromEnd,
			Name:              currNearestBed.Name,
			Score:             0,
			FieldsInitialized: 5}
		if currNearestBed.Score < contacts[j].B.ChromStart { // currNearestBed.Score is the abolute position of the TSS
			currAnswerB.Score = contacts[j].B.ChromStart - currNearestBed.Score
		} else {
			currAnswerB.Score = currNearestBed.Score - contacts[j].B.ChromStart
		}

		if currAnswerA.Score < currAnswerB.Score {
			currAnswerB.Score = currAnswerA.Score
			currAnswerB.Name = currAnswerA.Name
		} else if currAnswerB.Score < currAnswerA.Score {
			currAnswerA.Score = currAnswerB.Score
			currAnswerA.Name = currAnswerB.Name
		} //if scores are equal to each other they both retain their original closest gene assignment
		answer = append(answer, currAnswerA)
		answer = append(answer, currAnswerB)
	}

	return FillSpaceHiddenValue(answer, sizes)
}

// makeChromSlice is a helper function of Fill3dSpace that makes a list of all the chroms in the set of genes
func makeChromSlice(records []bed.Bed) []string {
	var chroms []string
	var inChroms bool

	for r := range records {
		if r == 0 {
			chroms = append(chroms, records[r].Chrom)
		} else {
			inChroms = false
			for n := range chroms {
				if records[r].Chrom == chroms[n] {
					inChroms = true
				}
			}
			if !inChroms {
				chroms = append(chroms, records[r].Chrom)
			}
		}
	}
	return chroms
}

// checkGeneFileForChrom is a helper function for fill3dSpace and checks the contacts bedpe to make sure that the chromosome is in the gene set as well
func checkGeneFileForChrom(a bedpe.BedPe, b []string) bool {
	chromInFile := false

	for r := range b {
		if b[r] == a.A.Chrom {
			chromInFile = true
		}
	}
	return chromInFile
}

// FillSpaceNoHiddenValue accepts a bed of single point tss positions and reference genome as a map[string]chromInfo.ChromInfo and returns a bed
// that assigns each genomic position to the nearest feature in the input bed using the Name field in the input bed.
// Absolute start positions of the original input bed entries are stored in the Score field of the output bed entries.
func FillSpaceNoHiddenValue(records []bed.Bed, genome map[string]chromInfo.ChromInfo) []bed.Bed {
	records = removeRecordsOnMissingChrom(records, genome)
	var answer = make([]bed.Bed, 0)
	if records == nil || len(records) == 0 {
		return records // return slice if it is empty or nil.
	}

	var currAnswer = bed.Bed{Chrom: records[0].Chrom, ChromStart: 0, ChromEnd: records[0].ChromEnd, Name: records[0].Name, Score: records[0].ChromStart, FieldsInitialized: 5}
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
func removeRecordsOnMissingChrom(records []bed.Bed, genome map[string]chromInfo.ChromInfo) []bed.Bed {
	var out = make([]bed.Bed, 0)
	var exist bool

	for i := range records {
		_, exist = genome[records[i].Chrom]
		if exist {
			out = append(out, records[i])
		}
	}
	return out
}

// FillSpaceHiddenValue accepts a slice of Bed structs and a reference genome as a map[string]chromInfo.ChromInfo and returns a
// slice of Bed structs that assigns each genomic position to the nearest feature in 3D space from the input bed, using
// the input bed scores to represent "hidden values", or the distance from that position to its nearest TSS in 3D space.
// bed that is passed here from Fill3dSpace has gene positions and bedpe positions already filled in
func FillSpaceHiddenValue(records []bed.Bed, genome map[string]chromInfo.ChromInfo) []bed.Bed {
	records = removeRecordsOnMissingChrom(records, genome)
	records = runUntilNoNewHidden(records)
	var violation = true
	for violation {
		records, violation = removeBedsWithNoTerritory(records)
	}
	var answer = make([]bed.Bed, 0)
	var threeDMidpoint int
	var currAnswer = bed.Bed{Chrom: records[0].Chrom, ChromStart: 0, ChromEnd: records[0].ChromEnd, Name: records[0].Name, FieldsInitialized: 4}

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

// runUntilNoNewHidden directs our mergeKeepLowNameAndScore function to keep updating until no changes occur,
// meaning everything is assigned the smallest possible distance to a gene
func runUntilNoNewHidden(records []bed.Bed) []bed.Bed {
	var newHidden bool
	records, newHidden = mergeKeepLowScoreAndName(records)
	if newHidden {
		runUntilNoNewHidden(records)
	}
	return records
}

// mergeKeepLowScoreAndName is a helper function for 3d great that takes in a slice of bed structs and determines if there's a mismatch in logic for hidden values.
// Hidden values should be the lowest posible distance between a bedpe point and a gene in 3d space. If there's a situation where
// the hidden value of a contact (A) is more than the distance to a neighbor (B) plus that neighbors hidden value (B's HV),
// then the hidden vale of A is recalculated as the distance from A to B plus B's hidden value.
func mergeKeepLowScoreAndName(records []bed.Bed) (out []bed.Bed, newHidden bool) {
	var dist int
	newHidden = false
	var outList = make([]bed.Bed, 0)
	bed.SortByCoord(records)
	var curr = records[0]
	for i := 1; i < len(records); i++ {
		if bed.Overlap(curr, records[i]) {
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

// helper fucntion of FillSpaceHiddenValue. Define the 3dMidpoint between bed entries A and B (with hidden values
// Ha and Hb) as (A.ChromEnd-Ha + B.ChromStart+Hb) / 2. If the midpoint is left of end of record A, then A is removed from the output
func removeBedsWithNoTerritory(records []bed.Bed) ([]bed.Bed, bool) {
	var answer []bed.Bed
	var threeDMidpoint int
	var violation = false
	for i := 1; i < len(records); i++ {
		if records[i-1].Chrom == records[i].Chrom {
			threeDMidpoint = (records[i-1].ChromEnd - records[i-1].Score + records[i].ChromStart + records[i].Score) / 2
			if threeDMidpoint < records[i-1].ChromEnd {
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
