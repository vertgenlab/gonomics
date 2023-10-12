package interval

import (
	"fmt"
	"github.com/vertgenlab/gonomics/numbers"
	"github.com/vertgenlab/gonomics/bed"
	"sort"
	"strings"
)

// CoordsToString takes in a interval and returns a string in the format of chr:start-end
func CoordsToString(i Interval) string {
	return fmt.Sprintf("%s:%v-%v", i.GetChrom(), i.GetChromStart(), i.GetChromEnd())
}

// IntervalSize calculates the size of the interval
func IntervalSize(i Interval) int {
	return i.GetChromEnd() - i.GetChromStart()
}

// OverlapSize calculates the size of the overlap betwen 2 intervals, assuming they are from the same genome
func OverlapSize(a, b Interval) int {
	if a.GetChrom() != b.GetChrom() {
		return 0
	}
	end := numbers.Min(a.GetChromEnd(), b.GetChromEnd())
	start := numbers.Max(a.GetChromStart(), b.GetChromStart())
	if end <= start {
		return 0
	} else {
		return end - start
	}

// IntervalSimilarity takes in two slices of interval and returns the proportion of elements in the first slice that overlap an element in the second slice, the proportion of elements in the second slice
// that overlap an element in the first slice, and the average of those two metrics (a metric of how similar the two interval sets are). Interval sets must not be self-overlapping
func IntervalSimilarity(a []Interval, b []Interval) (float64, float64, float64) {
	var overlapsA, overlapsB, allOverlapsA, allOverlapsB []Interval

	intervalMapA := BuildTree(a)
	intervalMapB := BuildTree(b)
	for _, bd := range b {
		overlapsA = Query(intervalMapA, bd, "any")
		for _, i := range overlapsA {
			allOverlapsA = append(allOverlapsA, i)
		}
	}
	for _, bd := range a {
		overlapsB = Query(intervalMapB, bd, "any")
		for _, i := range overlapsB {
			allOverlapsB = append(allOverlapsB, i)
		}
	}

	allUniqOverlapsA := Unique(allOverlapsA)
	allUniqOverlapsB := Unique(allOverlapsB)

	overlapPercA := float64(len(allUniqOverlapsA)) / float64(len(a))
	overlapPercB := float64(len(allUniqOverlapsB)) / float64(len(b))

	return overlapPercA, overlapPercB, (overlapPercA + overlapPercB) / 2
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

// compare is a helper function for SortByCoord. It returns -1 if a is before b, and 1 if b is before a. It will return 0 if the intervals are identical.haqerBedComps_test.go
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

// Unique takes in a slice of Interval and returns a slice of Interval with only unique interval entries (determined by coordinates). The first instance of each interval is kept.
// The output will be sorted by position
func Unique(regions []Interval) []Interval {
	var uniqueRegions []Interval

	SortByCoord(regions)
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
