package lift

import (
	"sort"
	"strings"
)

//allAreEqual returns true if two input slices of Lifts contain Lift entries that all return true for equal.
func allAreEqual(a []Lift, b []Lift) bool {
	if len(a) != len(b) {
		return false
	}
	for i := 0; i < len(a); i++ {
		if !equal(a[i], b[i]) {
			return false
		}
	}
	return true
}

//equal returns true if two input Lift entries have the same Chrom, ChromStart, and ChromEnd. False otherwise.
func equal(a Lift, b Lift) bool {
	if strings.Compare(a.GetChrom(), b.GetChrom()) != 0 {
		return false
	}
	if a.GetChromStart() != b.GetChromStart() {
		return false
	}
	if a.GetChromEnd() != b.GetChromEnd() {
		return false
	}
	return true
}

//sortBySize sorts Lift interfaces based on their size.
func sortBySize(l []Lift) {
	sort.Slice(l, func(i, j int) bool {return compareSize(l[i], l[j]) == -1})
}

//compareByCoord is the comparison for Lift interfaces based on genomic coordinate. Used in SortByCoord.
func compareByCoord(a Lift, b Lift) int {
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

//compareSize returns zero for Lifts with an equal length (ChromEnd - ChromStart) and otherwise returns the ordering of the two Lift structs.
func compareSize(a Lift, b Lift) int {
	sizeA := a.GetChromEnd() - a.GetChromStart()
	sizeB := b.GetChromEnd() - b.GetChromStart()
	if sizeA < sizeB {
		return -1
	}
	if sizeA > sizeB {
		return 1
	}
	return 0
}

//compareChromEnd returns zero for Lifts with an equal ChromEnd position and otherwise returns the ordering of the two Lift entries by chromEnd.
func compareChromEnd(a Lift, b Lift) int {
	if a.GetChromEnd() < b.GetChromEnd() {
		return -1
	}
	if a.GetChromEnd() > b.GetChromEnd() {
		return 1
	}
	return 0
}

//compareChromEndByChrom compares Lifts by chromosome and then by chromEnd, but not by chromStart.
func compareChromEndByChrom(a Lift, b Lift) int {
	chromComp := strings.Compare(a.GetChrom(), b.GetChrom())
	if chromComp != 0 {
		return chromComp
	}
	return compareChromEnd(a, b)
}

//SortByCoord sorts a slice of Lift interfaces by their coordinates.
func SortByCoord(b []Lift) {
	sort.Slice(b, func(i, j int) bool {return compareByCoord(b[i], b[j]) == -1})
}
