package interval

// Adapted from Mao, Eran & Luo 2019
// DOI: 10.1038/s41598-019-41451-3

import (
	"github.com/vertgenlab/gonomics/fileio"
	"sort"
)

type Interval interface {
	GetChrom() string
	GetChromStart() int
	GetChromEnd() int
	WriteToFileHandle(*fileio.EasyWriter)
	SetExclude()
}

type IntervalNode struct {
	val    Interval // only stored in leaf nodes
	data   []Interval
	xMid   int
	iLeft  []int
	iRight []int
	lChild *IntervalNode
	rChild *IntervalNode
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

	for i := 0; i < len(large); {
		if smallIdx > len(small)-1 {
			answer[i] = -1
			i++
			continue
		}
		if yLess(large[i], small[smallIdx]) || yEqual(large[i], small[smallIdx]) {
			answer[i] = smallIdx
			i++
		} else {
			smallIdx++
		}
	}

	return answer
}

func splitIntervalsByChr(intervals []Interval) map[string][]Interval {
	answer := make(map[string][]Interval)
	for i := 0; i < len(intervals); i++ {
		answer[intervals[i].GetChrom()] = append(answer[intervals[i].GetChrom()], intervals[i])
	}
	return answer
}

func BuildTree(intervals []Interval) map[string]*IntervalNode {
	answer := make(map[string]*IntervalNode)
	chrMap := splitIntervalsByChr(intervals)

	for idx, val := range chrMap {
		sortIntervals(val, yLess)
		answer[idx] = buildTree(val)
	}

	return answer
}

func buildTree(intervals []Interval) *IntervalNode {
	p := make([]Interval, len(intervals))
	pY := make([]Interval, len(intervals))
	copy(p, intervals)
	copy(pY, intervals)

	//TODO: Optimize the sorting to get rid of xLess sort up front by copying pLeft and pRight and feeding the original xLess sorted child into the recursive call
	// 1. Sort P by y-value, return an array of intervals Py.
	sortIntervals(p, xLess)
	//sortIntervals(pY, yLess)

	var answer *IntervalNode

	// 2. if P contains only one interval i then
	// 3. Creating a leaf node vleaf storing this interval. i.e., vleaf. interval = i
	if len(p) == 1 {
		answer = &IntervalNode{val: p[0], data: pY}
	} else {
		// 4. else
		// 5. Split P into Pleft and Pright, the subsets ≤ and > the median x-value xmid of P.
		var midPos int
		if len(p)%2 == 0 {
			midPos = (len(p) / 2) - 1
		} else {
			midPos = (len(p) - 1) / 2
		}

		pLeft := make([]Interval, len(p[:midPos+1]))
		pRight := make([]Interval, len(p[midPos+1:]))
		copy(pLeft, p[:midPos+1])
		copy(pRight, p[midPos+1:])

		// 6. Sort Pleft and Pright by y-value.
		sortIntervals(pLeft, yLess)
		sortIntervals(pRight, yLess)

		// 7. Create an FC-index Ileft from Py to Pleft.
		// 8. Create an FC-index Iright from Py to Pright
		iLeft := createFCIndex(pY, pLeft)
		iRight := createFCIndex(pY, pRight)

		// 9. vleft ← BuildRTFC(Pleft)
		// 10. vright ← BuildRTFC(Pright)
		vLeft := buildTree(pLeft)
		vRight := buildTree(pRight)

		// 11. Create a node v storing xmid, Ileft and Iright. v.x = xmid; v.data = Py;
		// v. lfc = Ileft; v. rfc = Iright; v. lchild = vleft; v. rchild = vright
		answer = &IntervalNode{
			data:   pY,
			xMid:   p[midPos].GetChromStart(),
			iLeft:  iLeft,
			iRight: iRight,
			lChild: vLeft,
			rChild: vRight,
		}
	}

	// 12. end if
	// 13. return v
	return answer
}

func Query(treeMap map[string]*IntervalNode, q Interval, relationship string) []Interval {
	var answer []Interval
	switch relationship {
	case "any":
		answer = append(answer, query(treeMap[q.GetChrom()], q, "o")...)
		answer = append(answer, query(treeMap[q.GetChrom()], q, "oi")...)
		answer = append(answer, query(treeMap[q.GetChrom()], q, "d")...)
		answer = append(answer, query(treeMap[q.GetChrom()], q, "di")...)
		answer = append(answer, query(treeMap[q.GetChrom()], q, "m")...)
		answer = append(answer, query(treeMap[q.GetChrom()], q, "mi")...)
		answer = append(answer, query(treeMap[q.GetChrom()], q, "s")...)
		answer = append(answer, query(treeMap[q.GetChrom()], q, "si")...)
		answer = append(answer, query(treeMap[q.GetChrom()], q, "f")...)
		answer = append(answer, query(treeMap[q.GetChrom()], q, "fi")...)
		answer = append(answer, query(treeMap[q.GetChrom()], q, "e")...)
	case "within":
		answer = append(answer, query(treeMap[q.GetChrom()], q, "d")...)
		answer = append(answer, query(treeMap[q.GetChrom()], q, "s")...)
		answer = append(answer, query(treeMap[q.GetChrom()], q, "f")...)
		answer = append(answer, query(treeMap[q.GetChrom()], q, "e")...)
	case "start":
		answer = append(answer, query(treeMap[q.GetChrom()], q, "s")...)
		answer = append(answer, query(treeMap[q.GetChrom()], q, "si")...)
		answer = append(answer, query(treeMap[q.GetChrom()], q, "e")...)
	case "end":
		answer = append(answer, query(treeMap[q.GetChrom()], q, "f")...)
		answer = append(answer, query(treeMap[q.GetChrom()], q, "fi")...)
		answer = append(answer, query(treeMap[q.GetChrom()], q, "e")...)
	case "equal":
		answer = query(treeMap[q.GetChrom()], q, "e")
	default:
		answer = query(treeMap[q.GetChrom()], q, relationship)
	}
	return answer
}

func query(tree *IntervalNode, q Interval, relationship string) []Interval {
	// 1. Transform interval query with respect to interval I and relationship R
	// to range query with x range [x1, x2] and y range [y1, y2] according to Table 1.
	x1, x2, y1, y2 := transform(q, relationship)

	var answer []Interval // 2. Initialize the result set S = {}

	// 3. Find the split node vsplit in range tree T where the paths
	// to x1 and x2 split, or the leaf where both paths end.
	vSplit := findSplit(x1, x2, tree)

	if vSplit == nil {
		return nil
	}

	if vSplit.val != nil { // 4. if vsplit is a leaf node then
		// 5. if vsplit. interval. x ∈ [x1, x2] and vsplit. interval. y ∈ [y1, y2] then
		if withinRange(tree.val, relationship, x1, x2, y1, y2) {
			// 6. Report the interval in vsplit, S = S ∪ {vsplit. interval}
			switch z := tree.val.(type) {
			case *AggregateInterval:
				answer = append(answer, query(z.tree[q.GetChrom()], q, relationship)...)
			default:
				answer = append(answer, z)
			}
		} // 7. end if
		return answer // 8. return S
	}

	// 9. Perform binary search on vsplit.data for y1 by y-value,
	// find the index i of the smallest element no less than y1 in vsplit.data.
	i := sort.Search(len(vSplit.data), func(i int) bool {
		return float32(vSplit.data[i].GetChromEnd()) >= y1
	})

	ri := i // save this value for later

	if i >= len(vSplit.data) {
		return nil
	}

	var v *IntervalNode
	v, i = vSplit.lChild, vSplit.iLeft[i] // 10. vsplit = vsplit.lchild, i = vsplit.lfc[i]

	for v.val == nil && i != -1 { // 11. while v is not a leaf and i ≠ 1 do // ADDED MODIFICATION TO 11
		if x1 <= float32(v.xMid) { // 12. if x1 ≤ v.x then
			j := v.iRight[i] // 13. j = v.rfc[i]

			// 14. while j ≠ −1 and v.rchild.data[j].y ≤ y2 and j ≤ v.rchild.data.size do
			for j != -1 && j < len(v.rChild.data) && float32(v.rChild.data[j].GetChromEnd()) <= y2 {
				if relationship == "m" || relationship == "mi" {
					if v.rChild.data[j].GetChromStart() != v.rChild.data[j].GetChromEnd() {
						switch z := v.rChild.data[j].(type) {
						case *AggregateInterval:
							answer = append(answer, query(z.tree[q.GetChrom()], q, relationship)...)
						default:
							answer = append(answer, z)
						}
					}
				} else {
					// 15. Report interval, S = S ∪ {v. rchild. data[j]}
					switch z := v.rChild.data[j].(type) {
					case *AggregateInterval:
						answer = append(answer, query(z.tree[q.GetChrom()], q, relationship)...)
					default:
						answer = append(answer, z)
					}
				}
				j++ // 16. j=j+1
			} // 17. end while

			i, v = v.iLeft[i], v.lChild // 18. i = v.lfc[i], v = v.lchild

		} else { // 19. else

			i, v = v.iRight[i], v.rChild //20. i=v.rfc[i], v=v.rchild

		} // 21. end if
	} // 22. end while

	// 23. if v is a leaf and v.interval.x ∈ [x1, x2] and v.interval.y ∈ [y1, y2] then
	if v.val != nil && withinRange(v.val, relationship, x1, x2, y1, y2) {
		answer = append(answer, v.val) // 24. Report the interval in v, S = S ∪ {v. interval}
	} // 25. end if

	v, i = vSplit.rChild, vSplit.iRight[ri] // 26. v = vsplit. rchild, i = vsplit. rcf[i]

	for v.val == nil && i != -1 { // 27. while v is not a leaf and i≠−1 do
		if x2 >= float32(v.xMid) { // 28. if x2 ≥ v. x then
			j := v.iLeft[i] // 29. j = v. lfc[i]

			// 30. while j ≠ −1 and v. lchild. data[j]. y ≤ y2 and j ≤ v. lchild. data. size do
			for j != -1 && j < len(v.lChild.data) && float32(v.lChild.data[j].GetChromEnd()) <= y2 {
				if relationship == "m" || relationship == "mi" {
					if v.lChild.data[j].GetChromStart() != v.lChild.data[j].GetChromEnd() {
						switch z := v.lChild.data[j].(type) {
						case *AggregateInterval:
							answer = append(answer, query(z.tree[q.GetChrom()], q, relationship)...)
						default:
							answer = append(answer, z)
						}
					}
				} else {
					switch z := v.lChild.data[j].(type) {
					case *AggregateInterval:
						answer = append(answer, query(z.tree[q.GetChrom()], q, relationship)...)
					default:
						answer = append(answer, z)
					}
				} // 31. Report interval, S = S ∪ {v. lchild. data[j]}
				j++ // 32. j=j+1
			} // 33. end while

			i, v = v.iRight[i], v.rChild // 34. i=v. rfc[i], v=v. rchild

		} else { // 35. else
			i, v = v.iLeft[i], v.lChild // 36. i = v. lfc[i], v = v. lchild
		} // 37. end if
	} // 38. end while

	// 39. if v is a leaf and v. interval. x ∈ [x1, x2] and v. interval. y ∈ [y1, y2] then
	if v.val != nil && withinRange(v.val, relationship, x1, x2, y1, y2) {
		switch z := v.val.(type) {
		case *AggregateInterval:
			answer = append(answer, query(z.tree[q.GetChrom()], q, relationship)...)
		default:
			answer = append(answer, z)
		} // 40. Report the interval in v, S = S ∪ {v. interval}
	} // 41. end if
	return answer // 42. return S
}

func withinRange(q Interval, relationship string, x1, x2, y1, y2 float32) bool {
	q1 := float32(q.GetChromStart())
	q2 := float32(q.GetChromEnd())
	if relationship == "m" || relationship == "mi" {
		if q1 == q2 {
			return false
		}
	}
	return (q1 >= x1 && q1 <= x2) && (q2 >= y1 && q2 <= y2)
}

func findSplit(x1, x2 float32, node *IntervalNode) *IntervalNode {

	for node.val == nil {
		if float32(node.xMid) < x1 {
			node = node.rChild
		} else if x2 < float32(node.xMid) {
			node = node.lChild
		} else {
			return node
		}
	}
	return nil
}
