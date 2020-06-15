package interval

import (
	"sort"
)

type Interval interface {
	GetChrom() string
	GetChromStart() int
	GetChromEnd() int
}

type IntervalNode struct {
	Interval
}

type IntervalSlice []Interval

func sortIntervals(s IntervalSlice, less func(a, b Interval) bool) {
	sort.Slice(s, func(i, j int) bool { return less(s[i], s[j]) })
}

func xLess(a, b Interval) bool {
	return a.GetChromStart() < b.GetChromStart()
}

func yLess(a, b Interval) bool {
	return a.GetChromEnd() < b.GetChromEnd()
}

func createFCIndex(large IntervalSlice, small IntervalSlice) []int {
	answer := make([]int, len(large))
	var largeIdx, smallIdx int

	for largeIdx = 0; largeIdx < len(large); largeIdx++ {

	}
	return answer
}

func BuildRTFC(intervals IntervalSlice) *IntervalNode {

	P := make(IntervalSlice, len(intervals))
	Py := make(IntervalSlice, len(intervals))
	copy(P, intervals)
	copy(Py, intervals)

	sortIntervals(P, xLess)
	sortIntervals(P, yLess)

	var answer, vLeaf *IntervalNode

	if len(P) == 1 {
		vLeaf = &IntervalNode{Interval: P[0]}
	} else {
		var midPos int
		if len(P)%2 == 0 {
			midPos = (len(P) / 2) - 1
		} else {
			midPos = (len(P) - 1) / 2
		}

		Pleft := P[:midPos+1]
		Pright := P[midPos+1:]

		sortIntervals(Pleft, yLess)
		sortIntervals(Pright, yLess)
	}

	return answer
}
