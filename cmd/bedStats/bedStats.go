//Command Group: BED Tools

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

func bedStats(s Settings) {
	var freq float64
	var matchedBeds, nonMatchedBeds []bed.Bed
	trueContacts := bedpe.Read(s.InContacts)
	testBeds := bed.Read(s.InTestBed)
	freq, matchedBeds, nonMatchedBeds = GeneAssignmentCheck(trueContacts, testBeds)

	bed.Write(s.OutMatched, matchedBeds)
	bed.Write(s.OutNonMatched, nonMatchedBeds)
	fmt.Println(freq)
}

// GeneAssignmentCheck returns information about whether gene assignments (based on coordinates and name fields) match
// entered true data set and how many regions.
// This can be used to compare any bedpe to any bed based off of overlap and name matching categories. Only the A foot of a bedpe will be checked if checkBothFeet is not set to true
func GeneAssignmentCheck(truth []bedpe.BedPe, test []bed.Bed) (regionMatchFrequency float64, matchesWithDistance []bed.Bed, nonMatches []bed.Bed) {
	var matches, truthAsBeds, nonMatchBeds []bed.Bed
	var matchCount, name, overlaps int
	var matchCountFreq float64
	var truthIntervals, currNearest []interval.Interval
	var trueBed, currNearestBed, matchedBed bed.Bed
	var names, chromList []string
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

	for t := range mergedTruthBeds {
		if t == 0 {
			chromList = append(chromList, mergedTruthBeds[t].Chrom)
		} else if t != 0 && mergedTruthBeds[t].Chrom != mergedTruthBeds[t-1].Chrom {
			chromList = append(chromList, mergedTruthBeds[t].Chrom)
		}
		truthIntervals = append(truthIntervals, mergedTruthBeds[t])
	}

	log.Print(chromList)

	truthTree := interval.BuildTree(truthIntervals)

	for currTestBed := range test {
		if currTestBed > 0 && bed.Equal(test[currTestBed-1], test[currTestBed]) {
			continue
		}
		if currTestBed > 0 && len(currNearest) > 0 && !matched {
			nonMatchBeds = append(nonMatchBeds, test[currTestBed-1])
		}
		matched = false
		currNearest = interval.Query(truthTree, test[currTestBed], "any")
		if len(currNearest) == 0 { //we can have regions that don't have a contact in them, we will ignore those and not include them in our counts
			continue
		} else {
			overlaps++
		}
		for j = range currNearest {
			if matched {
				continue
			}
			currNearestBed = currNearest[j].(bed.Bed)
			names = strings.Split(currNearestBed.Name, ",")
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
				if names[name] == test[currTestBed].Name {
					matchCount++
					matched = true
					matchedBed = bed.Bed{Chrom: test[currTestBed].Chrom,
						ChromStart:        test[currTestBed].ChromStart,
						ChromEnd:          test[currTestBed].ChromEnd,
						Name:              test[currTestBed].Name + "," + names[name],
						FieldsInitialized: 7,
						Annotation:        append(matchedBed.Annotation, currNearestBed.Annotation[name])}
					matches = append(matches, matchedBed)
				} else {
					continue
				}
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
			"bedStats true.bedpe test.bed matched.bed nonMatched.bed \n" +
			"bedStats - takes a bedpe containing true contacts from empirical data which assigns regions of the genome \n" +
			"to putative target genes and compares an output from assignGenomeSpace command to determine how \n" +
			"accurate the nearest gene in 3d space calculation was. Writes out a frequency of correct assignments \n" +
			"and outputs a bed containing the matching regions and regions that overlapped the true data set but didn't match.\n" +
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
	bedStats(s)
}
