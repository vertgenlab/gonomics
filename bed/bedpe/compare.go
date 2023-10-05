package bedpe

import (
	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/interval"
	"github.com/vertgenlab/gonomics/numbers"
	"log"
	"strings"
)

func AllAreEqual(a []BedPe, b []BedPe) bool {
	if len(a) != len(b) {
		return false
	}
	for i := range a {
		if !Equal(a[i], b[i]) {
			return false
		}
	}
	return true
}

// Equal returns true if two input BedPe entries have the same coordinates for both regions. False otherwise.
func Equal(a BedPe, b BedPe) bool {
	if !bed.Equal(a.A, b.A) {
		return false
	}
	if !bed.Equal(a.B, b.B) {
		return false
	}
	return true
}

// GeneAssignmentCheck returns information about whether gene assignments (based on coordinates and name fields) match
// entered true data set and how many regions.
// This can be used to compare any bedpe to any bed based off of overlap and name matching categories. Only the A foot of a bedpe will be checked if checkBothFeet is not set to true
func GeneAssignmentCheck(truth []BedPe, test []bed.Bed) (regionMatchFrequency float64, matchesWithDistance []bed.Bed) {
	var matches, truthAsBeds []bed.Bed
	var matchCount, nonMatchCount, name int
	var matchCountFreq float64
	var truthIntervals, currNearest []interval.Interval
	var trueBed, currNearestBed, matchedBed bed.Bed
	var names []string
	var matched bool
	var chromList []string
	var j int

	annotateTruthFeetDist(truth)

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
		} else if mergedTruthBeds[t].Chrom != mergedTruthBeds[t-1].Chrom {
			chromList = append(chromList, mergedTruthBeds[t].Chrom)
		}
		truthIntervals = append(truthIntervals, mergedTruthBeds[t])
	}

	truthTree := interval.BuildTree(truthIntervals)

	for currTestBed := range test {
		matched = false
		currNearest = interval.Query(truthTree, test[currTestBed], "any")
		if len(currNearest) == 0 { //we can have regions that don't have a contact in them, we will ignore those and not include them in our counts
			continue
		}

		for j = range currNearest {
			currNearestBed = currNearest[j].(bed.Bed)
			names = strings.Split(currNearestBed.Name, ",")
			for name = range names {
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
						Name:              test[currTestBed].Name,
						FieldsInitialized: 7,
						Annotation:        append(matchedBed.Annotation, currNearestBed.Annotation[name])}
					matches = append(matches, matchedBed)
				} else {
					continue
				}
			}
		}
		if !matched {
			nonMatchCount++
		}
	}

	log.Printf("matches: %v, nonMatches: %v", matchCount, nonMatchCount)
	matchCountFreq = float64(matchCount) / float64(nonMatchCount+matchCount)
	return matchCountFreq, matches
}

func annotateTruthFeetDist(b []BedPe) {
	var dist int
	for i := range b {
		b[i].A.FieldsInitialized = 11
		dist = numbers.Max(b[i].A.ChromStart, b[i].B.ChromStart) - numbers.Min(b[i].A.ChromStart, b[i].B.ChromStart)
		b[i].A.Annotation = append(b[i].A.Annotation, fileio.IntToString(dist))
	}
}
