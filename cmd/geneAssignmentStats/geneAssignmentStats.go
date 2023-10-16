//Command Group: BED Tools

// geneAssignmentStats compares a bedpe containing a gene symbol in the name field to a test set of output from
// assignGenomeSpace command which assigns every base in the genome a closest gene either from proximity or from a 3d contact map.
package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/bed/bedpe"
	"github.com/vertgenlab/gonomics/interval"
	"log"
	"strings"
)

type Settings struct {
	InContacts    string
	InTestBed     string
	OutMatched    string
	OutNonMatched string
}

func geneAssignmentStats(s Settings) {
	var freq float64
	var matchedBeds, nonMatchedBeds []bed.Bed
	trueContacts := bedpe.Read(s.InContacts)
	testBeds := bed.Read(s.InTestBed)
	freq, matchedBeds, nonMatchedBeds = GeneAssignmentCheckGuidePers(trueContacts, testBeds)

	bed.Write(s.OutMatched, matchedBeds)
	bed.Write(s.OutNonMatched, nonMatchedBeds)
	fmt.Println(freq)
}

//// GeneAssignmentCheck returns information about whether gene assignments (based on coordinates and name fields) match
//// entered true data set and how many regions.
//// This can be used to compare any bedpe to any bed based off of overlap and name matching categories. Only the A foot of a bedpe will be checked if checkBothFeet is not set to true
//func GeneAssignmentCheck(truth []bedpe.BedPe, test []bed.Bed) (regionMatchFrequency float64, matchesWithDistance []bed.Bed, nonMatches []bed.Bed) {
//	var matches, truthAsBeds, nonMatchBeds []bed.Bed
//	var matchCount, name, overlaps int
//	var matchCountFreq float64
//	var truthIntervals, currNearest []interval.Interval
//	var trueBed, currNearestBed, matchedBed bed.Bed
//	var names []string
//	var matched bool
//	var j int
//
//	bedpe.AnnotateFeetDist(truth)
//
//	for currTruth := range truth {
//		trueBed = bed.Bed{
//			Chrom:             truth[currTruth].A.Chrom,
//			ChromStart:        truth[currTruth].A.ChromStart,
//			ChromEnd:          truth[currTruth].A.ChromEnd,
//			Name:              truth[currTruth].A.Name,
//			Annotation:        truth[currTruth].A.Annotation,
//			FieldsInitialized: 7}
//		truthAsBeds = append(truthAsBeds, trueBed)
//	}
//
//	mergedTruthBeds := bed.MergeBedsKeepNamesAndAnnotations(truthAsBeds)
//
//	bed.Write("mergedBeds.GasperiniContacts.bed", mergedTruthBeds)
//
//	for t := range mergedTruthBeds {
//		truthIntervals = append(truthIntervals, mergedTruthBeds[t])
//	}
//
//	log.Print(len(truthIntervals))
//
//	truthTree := interval.BuildTree(truthIntervals)
//
//	log.Print(len(test))
//	for currTestBed := range test {
//		matched = false
//		currNearest = interval.Query(truthTree, test[currTestBed], "any")
//		if len(currNearest) == 0 { //we can have regions that don't have a contact in them, we will ignore those and not include them in our counts
//			continue
//		} else {
//			overlaps += len(currNearest)
//		}
//		for j = range currNearest {
//			if matched {
//				continue
//			}
//			currNearestBed = currNearest[j].(bed.Bed)
//			names = strings.Split(currNearestBed.Name, ",")
//			for name = range names {
//				if matched {
//					continue
//				}
//				matchedBed = bed.Bed{Chrom: "",
//					ChromStart:        0,
//					ChromEnd:          0,
//					Name:              "",
//					FieldsInitialized: 7,
//					Annotation:        []string{}}
//				if names[name] == test[currTestBed].Name {
//					matchCount++
//					matched = true
//					matchedBed = bed.Bed{Chrom: test[currTestBed].Chrom,
//						ChromStart:        test[currTestBed].ChromStart,
//						ChromEnd:          test[currTestBed].ChromEnd,
//						Name:              test[currTestBed].Name + "," + names[name],
//						FieldsInitialized: 7,
//						Annotation:        append(matchedBed.Annotation, currNearestBed.Annotation[name])}
//					matches = append(matches, matchedBed)
//				} else {
//					nonMatchBeds = append(nonMatchBeds, test[currTestBed])
//					continue
//				}
//			}
//		}
//	}
//	//divided the number of regions with overlap and matching names by the number of regions in the true data set
//	matchCountFreq = float64(matchCount) / float64(len(mergedTruthBeds))
//	log.Printf("Matched: %v, Total: %v, Number of regions in test that overlapped true: %v", matchCount, len(mergedTruthBeds), overlaps)
//	return matchCountFreq, matches, nonMatchBeds
//}

// GeneAssignmentCheck returns information about whether gene assignments (based on coordinates and name fields) match
// entered true data set and how many regions.
// This can be used to compare any bedpe to any bed based off of overlap and name matching categories. Only the A foot of a bedpe will be checked if checkBothFeet is not set to true
func GeneAssignmentCheckGuidePers(truth []bedpe.BedPe, test []bed.Bed) (regionMatchFrequency float64, matchesWithDistance []bed.Bed, nonMatches []bed.Bed) {
	var matches, truthAsBeds, nonMatchBeds []bed.Bed
	var matchCount, name, overlaps int
	var matchCountFreq float64
	var testIntervals, currNearest []interval.Interval
	var trueBed, currNearestBed, matchedBed, nonMatchBed bed.Bed
	var names []string
	var matched bool
	var j int

	bedpe.AnnotateFeetDist(truth)

	for currTruth := range truth {
		trueBed = bed.Bed{
			Chrom:             truth[currTruth].A.Chrom,
			ChromStart:        truth[currTruth].A.ChromStart,
			ChromEnd:          truth[currTruth].A.ChromEnd,
			Name:              truth[currTruth].A.Name,
			Annotation:        truth[currTruth].A.Annotation,
			FieldsInitialized: 7}
		truthAsBeds = append(truthAsBeds, trueBed)
	}

	mergedTruthBeds := bed.MergeBedsKeepNamesAndAnnotations(truthAsBeds)

	for t := range test {
		testIntervals = append(testIntervals, test[t])
	}

	testTree := interval.BuildTree(testIntervals)

	for currTrue := range mergedTruthBeds {
		matched = false
		currNearest = interval.Query(testTree, mergedTruthBeds[currTrue], "any")
		if len(currNearest) == 0 {
			log.Fatal("No overlap found for guide")
		}
		for j = range currNearest {
			currNearestBed = currNearest[j].(bed.Bed)
			names = strings.Split(mergedTruthBeds[currTrue].Name, ",")
			for name = range names {
				if matched {
					continue
				}
				matchedBed = bed.Bed{Chrom: "",
					ChromStart:        0,
					ChromEnd:          0,
					Name:              "",
					FieldsInitialized: 7,
					Annotation:        []string{}}
				if names[name] == currNearestBed.Name {
					matchCount++
					matched = true
					matchedBed = bed.Bed{Chrom: mergedTruthBeds[currTrue].Chrom,
						ChromStart:        mergedTruthBeds[currTrue].ChromStart,
						ChromEnd:          mergedTruthBeds[currTrue].ChromEnd,
						Name:              names[name] + "," + currNearestBed.Name,
						FieldsInitialized: 7,
						Annotation:        append(matchedBed.Annotation, mergedTruthBeds[currTrue].Annotation[name])}
					matches = append(matches, matchedBed)
				}
			}
			if !matched {
				nonMatchBed = bed.Bed{Chrom: mergedTruthBeds[currTrue].Chrom,
					ChromStart:        mergedTruthBeds[currTrue].ChromStart,
					ChromEnd:          mergedTruthBeds[currTrue].ChromEnd,
					Name:              names[name] + "," + currNearestBed.Name,
					FieldsInitialized: 7,
					Annotation:        append(matchedBed.Annotation, mergedTruthBeds[currTrue].Annotation[name])}
				nonMatchBeds = append(nonMatchBeds, nonMatchBed)
			}
		}
	}
	//divided the number of regions with overlap and matching names by the number of regions in the true data set
	matchCountFreq = float64(matchCount) / float64(len(mergedTruthBeds))
	log.Printf("Matched: %v, Total: %v, Number of regions in test that overlapped true: %v", matchCount, len(mergedTruthBeds), overlaps)
	return matchCountFreq, matches, nonMatchBeds
}

func usage() {
	fmt.Print(
		"Usage: \n" +
			"geneAssignmentStats true.bedpe test.bed matched.bed nonMatched.bed \n" +
			"geneAssignmentStats - takes a bedpe containing true contacts from empirical data which assigns regions of " +
			"the genome (guide coordinates kept as the first record in the bedpe) \n" +
			"to putative target genes (coordinates kept as the second set of coordinates in the bedpe and the name of " +
			"the bedpe file will correspond to the gene symbol) and compares an output from assignGenomeSpace command to determine how \n" +
			"accurate the nearest gene assignment is. Writes out a frequency of correct assignments \n" +
			"and outputs a bed containing the matching regions where the name field should have the matching gene symbols " +
			"in comma separated format and a bed file of guide coordinates which have a name field of their gene " +
			"assignment and fthe gene assignment that was incorrect in the input test bed file.\n" +
			"options:\n")
	flag.PrintDefaults()
}
func main() {
	var expectedNumArgs int = 4

	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime | log.Lshortfile)
	flag.Parse()
	if len(flag.Args()) != expectedNumArgs {
		flag.Usage()
		log.Fatalf("Error: expecting %d arguments, but got %d\n",
			expectedNumArgs, len(flag.Args()))
	}
	inContacts := flag.Arg(0)
	inTest := flag.Arg(1)
	matchedBeds := flag.Arg(2)
	nonMatchedBeds := flag.Arg(3)

	s := Settings{
		InContacts:    inContacts,
		InTestBed:     inTest,
		OutMatched:    matchedBeds,
		OutNonMatched: nonMatchedBeds,
	}
	geneAssignmentStats(s)
}
