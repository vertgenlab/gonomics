package interval

import (
	"fmt"
	"github.com/vertgenlab/gonomics/bed"
	"sort"
	"strings"
)

// CoordsToString takes in a interval and returns a string in the format of chr:start-end
func CoordsToString(i Interval) string {
	return fmt.Sprintf("%s:%v-%v", i.GetChrom(), i.GetChromStart(), i.GetChromEnd())
}

// IntervalSimilarity takes in two slices of interval and returns the proportion of elements in the first slice that overlap an element in the second slice, the proportion of elements in the second slice
// that overlap an element in the first slice, and the average of those two metrics (a metric of how similar the two interval sets are). Interval sets must not be self-overlapping
func IntervalSimilarity(a []Interval, b []Interval) (float64, float64, float64) {
	var overlapSmall, overlapLarge, allLargeOverlaps, allSmallOverlaps []Interval

	smallIntervalMap := BuildTree(a)
	largeIntervalMap := BuildTree(b)
	for _, bd := range b {
		overlapSmall = Query(smallIntervalMap, bd, "any")
		for _, i := range overlapSmall {
			allSmallOverlaps = append(allSmallOverlaps, i)
		}
	}
	for _, bd := range a {
		overlapLarge = Query(largeIntervalMap, bd, "any")
		for _, i := range overlapLarge {
			allLargeOverlaps = append(allLargeOverlaps, i)
		}
	}

	SortByCoord(allSmallOverlaps)
	SortByCoord(allLargeOverlaps)

	allSmallOverlapsUniq := Unique(allSmallOverlaps)
	allLargeOverlapsUniq := Unique(allLargeOverlaps)

	smallOverlapPerc := float64(len(allSmallOverlapsUniq)) / float64(len(a))
	bigOverlapPerc := float64(len(allLargeOverlapsUniq)) / float64(len(b))

	return smallOverlapPerc, bigOverlapPerc, (smallOverlapPerc + bigOverlapPerc) / 2
}

// BedSliceToIntervalMap takes in a slice of bed and returns a map of string -- *IntervalNode which can be used for interval.Query
func BedSliceToIntervalMap(inBed []bed.Bed) map[string]*IntervalNode {
	var smallWindowIntervals []Interval

	for bd := range inBed {
		smallWindowIntervals = append(smallWindowIntervals, inBed[bd])
	}
	return BuildTree(smallWindowIntervals)
}

// AreEqual takes in two intervals and returns a bool of whether they have the same coordinates or not
func AreEqual(a Interval, b Interval) bool {
	if a.GetChrom() == b.GetChrom() && a.GetChromStart() == b.GetChromStart() && a.GetChromEnd() == b.GetChromEnd() {
		return true
	} else {
		return false
	}
}

// SortByCoord takes in a slice of Interval and sorts the slice according to coordinates
func SortByCoord(in []Interval) {
	sort.Slice(in, func(i, j int) bool { return compare(in[i], in[j]) == -1 })
}

func compare(a Interval, b Interval) int {
	chromComp := strings.Compare(a.GetChrom(), b.GetChrom())
	if chromComp != 0 {
		return chromComp
	}
	if a.GetChromStart() < b.GetChromStart() {
		return -1
	}
	if a.GetChromStart() > b.GetChromStart() {
		return 1
	}
	if a.GetChromEnd() < b.GetChromEnd() {
		return -1
	}
	if a.GetChromEnd() > b.GetChromEnd() {
		return 1
	}
	return 0
}

// Unique takes in a slice of Interval and returns a slice of Interval with only unique interval entries (determined by coordinates). The last instance of each interval is kept.
func Unique(regions []Interval) []Interval {
	var uniqueRegions []Interval

	currRegion := regions[0]
	for _, i := range regions {
		if AreEqual(currRegion, i) {
			continue
		}
		uniqueRegions = append(uniqueRegions, currRegion)
		currRegion = i
	}
	uniqueRegions = append(uniqueRegions, currRegion)
	return uniqueRegions
}

// BedSliceToIntervals takes in a slice of beds and returns a slice of intervals. Useful before using BuildTree
func BedSliceToIntervals(inBed []bed.Bed) []Interval {
	var intervalSlice []Interval

	for bd := range inBed {
		intervalSlice = append(intervalSlice, inBed[bd])
	}
	return intervalSlice
}
