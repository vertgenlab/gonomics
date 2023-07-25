package bedpe

import (
	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/chromInfo"
	"github.com/vertgenlab/gonomics/interval"
	"log"
	//TODO: Reintroduce when bug is fixed "log"
	//"fmt"
)

// Fill3dSpace takes in a set of Bedpe format contact points, single-base pair []bed.Bed representing
// trancription start sites, and a reference genome as a map[string]chromInfo.ChromInfo. This program
// returns a []bed.Bed representing the closest TSS in 3D space for each position in the reference genome,
// using the Bedpe contacts to determine the closest TSS in 3D space.
// Note that the input TSS must have 0 in score fields.
func Fill3dSpace(contacts []BedPe, tss []bed.Bed, sizes map[string]chromInfo.ChromInfo) []bed.Bed {
	var chromInFile bool
	geneChroms := makeChromSlice(tss) //makes a list of the chroms in the gene file
	answer := make([]bed.Bed, len(tss))
	copy(answer, tss)
	var currNearestBed, currAnswerA, currAnswerB bed.Bed
	var closest1dGenesIntervals, currNearest []interval.Interval
	closest1dGene := bed.FillSpaceNoHiddenValue(tss, sizes)

	for i := range closest1dGene {
		closest1dGenesIntervals = append(closest1dGenesIntervals, closest1dGene[i])
	}
	closest1dGeneTree := interval.BuildTree(closest1dGenesIntervals)

	contactsToMidpoints(contacts)

	// this for loop finds the nearest gene and hidden value for each bedpe foot midpoint
	for j := range contacts {
		//this is just a check to make sure that for any new chromosomes we encounter that they exist in the gene file, so we don't throw an error.
		if j == 0 {
			chromInFile = checkGeneFileForChrom(contacts[j], geneChroms)
			if !chromInFile {
				continue
			}
		} else if contacts[j].A.Chrom != contacts[j-1].A.Chrom { // only check when we encounter a new chrom in the bedpe file
			chromInFile = checkGeneFileForChrom(contacts[j], geneChroms)
			if !chromInFile {
				continue
			}
		}
		currNearest = interval.Query(closest1dGeneTree, contacts[j].A, "any")
		if len(currNearest) > 1 || len(currNearest) == 0 {
			log.Fatal("Space Filled bed should return one nearest bed entry.")
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

	return bed.FillSpaceHiddenValue(answer, sizes)
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
func checkGeneFileForChrom(a BedPe, b []string) bool {
	chromInFile := false

	for r := range b {
		if b[r] == a.A.Chrom {
			chromInFile = true
		}
	}
	return chromInFile
}
