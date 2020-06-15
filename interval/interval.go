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

func Query(tree *IntervalNode, q Interval, relationship string) []Interval {
	var answer []Interval
	// TODO: build in logic for: any, within, start, end, equal
	// see table 7 for relationships
	switch relationship {
	case "any":
	case "within":
	case "start":
	case "end":
	case "equal":
	default:
		answer = query(tree, q, relationship)
	}
	return answer
}

func query(tree *IntervalNode, q Interval, relationship string) []Interval {
	var answer []Interval
	x1, x2, y1, y2 := transform(q, relationship)

	vSplit := findSplit()

	if vSplit == nil {
		return nil
	}

	if vSplit.val != nil {
		if withinRange(tree.val, x1, x2, y1, y2) {
			answer = append(answer, tree.val)
		}
	}
	// TODO: pickup on algorithm 2 line 9

	return answer
}

func withinRange(q Interval, x1, x2, y1, y2 float32) bool {
	q1 := float32(q.GetChromStart())
	q2 := float32(q.GetChromEnd())
	return (q1 > x1 && q1 < x2) && (q2 > y1 && q2 < y2)
}

func findSplit() *IntervalNode {

}