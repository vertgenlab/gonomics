package lift

import (
	"fmt"
	"github.com/vertgenlab/gonomics/numbers"
	"strings"
)

//OverlapCount returns the number of elements from list one that have any overlap with list two. Answers range from 0 to len(a).
//Input Lift slices must be presorted with SortByCoord.
func OverlapCount(a []Lift, b []Lift) int {
	var count int = 0
	var aIndex, bIndex int

	for aIndex < len(a) && bIndex < len(b) {
		if overlap(a[aIndex], b[bIndex]) {
			count++
			aIndex++
		} else if compareChromEndByChrom(a[aIndex], b[bIndex]) < 0 {
			aIndex++
		} else {
			bIndex++
		}
	}
	return count
}

// overlapProbability calculates the probability that an element of len 'length' overlaps a set of genomic 'elements' in a genome represented by 'noGapRegions'.
// tempElements and tempNoGap represent cloned slices of elements and noGapRegions which are used for memory optimization.
func overlapProbability(elements []Lift, tempElements []Lift, length int, noGapRegions []Lift, tempNoGap []Lift) float64 {
	subtractFromCoord(elements, length-1, 0, tempElements)
	subtractFromCoord(noGapRegions, 0, length-1, tempNoGap)
	return float64(overlapLengthSum(tempElements, tempNoGap)) / float64(totalSize(tempNoGap))
}

//overlapLength returns the number of bases for which two Lift entries overlap.
func overlapLength(a Lift, b Lift) int {
	if !overlap(a, b) {
		return 0
	}
	end := numbers.Min(a.GetChromEnd(), b.GetChromEnd())
	start := numbers.Max(a.GetChromStart(), b.GetChromStart())
	return end - start
}

//overlapLengthSum calculates the total number of overlapping bases between two sets of Lift elements.
//Input Lift slices must be presorted with sortByCoord
func overlapLengthSum(a []Lift, b []Lift) int {
	var sum int = 0
	var aIndex, bIndex, oLen int
	for aIndex < len(a) && bIndex < len(b) {
		oLen = overlapLength(a[aIndex], b[bIndex])
		if oLen != 0 {
			sum += oLen
		}
		if compareChromEndByChrom(a[aIndex], b[bIndex]) < 0 {
			aIndex++
		} else {
			bIndex++
		}
	}
	return sum
}

//totalSize returns the total length covered in a slice of Lift entries.
func totalSize(b []Lift) int {
	var ans, curLen int
	for i := 0; i < len(b); i++ {
		curLen = b[i].GetChromEnd() - b[i].GetChromStart()
		ans += curLen
	}
	return ans
}

// helperfunction that returns an int corresponding to the shortest distance between chromStart and chromEnd from an input set of Lift entries.
func findShortestLength(b []Lift) int {
	var minLength int = b[0].GetChromEnd() - b[0].GetChromStart()
	for i := 1; i < len(b); i++ {
		if b[i].GetChromEnd()-b[i].GetChromStart() < minLength {
			minLength = b[i].GetChromEnd() - b[i].GetChromStart()
		}
	}
	return minLength
}

// helperfunction that returns an int corresponding to the longest distance between chromStart and chromEnd from an input set of Lift entries.
func findLargestLength(b []Lift) int {
	var maxLength int = b[0].GetChromEnd() - b[0].GetChromStart()
	for i := 1; i < len(b); i++ {
		if b[i].GetChromEnd()-b[i].GetChromStart() > maxLength {
			maxLength = b[i].GetChromEnd() - b[i].GetChromStart()
		}
	}
	return maxLength
}

//IsSelfOverlapping returns true if any elements in a Lift slice overlap another element in the same slice. False otherwise.
//b must be presorted with SortByCoord. Verbose > 0 reveals debug prints.
func IsSelfOverlapping(b []Lift, verbose int) bool {
	for i := 0; i < len(b)-1; i++ {
		if overlap(b[i], b[i+1]) {
			if verbose > 0 {
				fmt.Printf("first Lift: %v\n", b[i])
				fmt.Printf("second Lift: %v\n", b[i+1])
			}
			return true
		}
	}
	return false
}

func overlap(alpha Lift, beta Lift) bool {
	if (numbers.Max(alpha.GetChromStart(), beta.GetChromStart()) < numbers.Min(alpha.GetChromEnd(), beta.GetChromEnd())) && strings.Compare(alpha.GetChrom(), beta.GetChrom()) == 0 {
		return true
	} else {
		return false
	}
}
