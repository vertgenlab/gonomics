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
	val 	Interval
	data 	[]Interval
	xMid 	int
	iLeft 	[]int
	iRight 	[]int
	lChild 	*IntervalNode
	rChild 	*IntervalNode
}

func sortIntervals(s []Interval, less func(a, b Interval) bool) {
	sort.Slice(s, func(i, j int) bool { return less(s[i], s[j]) })
}

func xLess(a, b Interval) bool {
	return a.GetChromStart() < b.GetChromStart()
}

func yLess(a, b Interval) bool {
	return a.GetChromEnd() < b.GetChromEnd()
}

func yEqual(a, b Interval) bool {
	return a.GetChromEnd() == b.GetChromEnd()
}

func createFCIndex(large []Interval, small []Interval) []int {
	answer := make([]int, len(large))
	var smallIdx int

	for idx, val := range large {
		if smallIdx > len(small) - 1 {
			answer[idx] = -1
			 continue
		}
		if yLess(val, small[smallIdx]) || yEqual(val, small[smallIdx]) {
			answer[idx] = smallIdx
		} else {
			smallIdx++
		}
	}

	return answer
}

func BuildTree(intervals []Interval) *IntervalNode {

	p := make([]Interval, len(intervals))
	pY := make([]Interval, len(intervals))
	copy(p, intervals)
	copy(pY, intervals)

	sortIntervals(p, xLess)
	sortIntervals(p, yLess)

	var answer *IntervalNode

	if len(p) == 1 {
		answer = &IntervalNode{val: p[0]}
	} else {
		var midPos int
		if len(p)%2 == 0 {
			midPos = (len(p) / 2) - 1
		} else {
			midPos = (len(p) - 1) / 2
		}

		pLeft := p[:midPos+1]
		pRight := p[midPos+1:]

		sortIntervals(pLeft, yLess)
		sortIntervals(pRight, yLess)

		iLeft := createFCIndex(pY, pLeft)
		iRight := createFCIndex(pY, pRight)

		vLeft := BuildTree(pLeft)
		vRight := BuildTree(pRight)

		answer = &IntervalNode{
			data: 	pY,
			xMid:	p[midPos].GetChromStart(),
			iLeft:	iLeft,
			iRight: iRight,
			lChild: vLeft,
			rChild: vRight,
		}
	}

	return answer
}

func Query(tree *IntervalNode, query Interval, relationship string) []Interval {

}