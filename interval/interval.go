package interval

import (
	"sort"
)

type Interval interface {
	GetChrom() string
	GetChromStart() int
	GetChromEnd() int
}

type IntervalSlice []Interval

func sortIntervals(s IntervalSlice) {
	less := func(i, j int) bool {
		return s[i].GetChromStart() < s[j].GetChromStart() ||
			(s[i].GetChromStart() == s[j].GetChromStart() && s[i].GetChromEnd() <= s[j].GetChromEnd())
	}
	sort.Slice(s, less)
}

// For each interval (l) in IntervalSlice
// holds all intervals (a) in IntervalSlice
// in which the len of a < len of l such that
// such that smallerIntervals is a map[l][]a
// Note: []a is sorted by length
type smallerIntervals map[Interval]IntervalSlice

// Generates smallerIntervals map for a input list of intervals
// Removes each element put in smallerIntervals from s
// Note: must input a sorted interval slice
func computeSmallerIntervals(s IntervalSlice) (IntervalSlice, smallerIntervals) {
	smaller := make(smallerIntervals)
	var keptIntervals int = len(s)

	for i := 0; i < len(s)-1; i++ {
		// If a two intervals are found to have the same start
		if s[i].GetChromStart() == s[i+1].GetChromStart() {
			var currSmaller IntervalSlice

			// Go through slice until the left endpoints do not equal
			// first check is to make sure we dont run off the slice
			for ; len(s) > i+1 && s[i].GetChromStart() == s[i+1].GetChromStart(); i++ {
				// Add s[i] to currSmaller and remove from s
				currSmaller = append(currSmaller, s[i])
				s[i] = nil
				keptIntervals--
			}
			// Flip currSmaller so that intervals are sorted by descending order
			reverse(currSmaller)
			// Make the entry in the map and link to currSmaller
			smaller[s[i]] = currSmaller
		}
	}

	// Recast s into new interval slice with the nil values removed
	newIntervals := make(IntervalSlice, keptIntervals)
	var currIdx int = 0

	for i := 0; i < len(s); i++ {
		if s[i] != nil {
			newIntervals[currIdx] = s[i]
			currIdx++
		}
	}

	return newIntervals, smaller
}

// Flips slice
func reverse(s IntervalSlice) {
	for i, j := 0, len(s)-1; i < j; i, j = i+1, j-1 {
		s[i], s[j] = s[j], s[i]
	}
}

// s must be sorted and have gone through computeSmallerIntervals
func ConstructTree(s IntervalSlice) []*IntervalNode {
	answer := make([]*IntervalNode, len(s)+1)
	dummy := &IntervalNode{nil, nil, nil, make([]*IntervalNode, 0)}
	answer[0] = dummy

	// Add first node to dummy
	answer[1] = &IntervalNode{Interval: s[0], Parent: dummy, Children: make([]*IntervalNode, 0)}
	dummy.Children = append(dummy.Children, answer[1])

	for i := 1; i < len(s); i++ {

	}
	return answer
}
