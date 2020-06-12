package interval

type IntervalNode struct {
	Interval
	Parent   *IntervalNode
	LeftSib  *IntervalNode
	Children []*IntervalNode
}

func traverseTree(v *IntervalNode, stack *[]Interval, q Interval, smaller smallerIntervals) {
	// Output current node
	*stack = append(*stack, v.Interval)

	// Check all nodes in smaller intervals of v for any that intersect the query
	var w Interval
	for _, w = range smaller[v.Interval] {
		if w.GetChromEnd() < q.GetChromEnd() {
			break
		}
		*stack = append(*stack, w)
	}

	// Calculate the rightmost path through the tree
	// Definition 1. The rightmost path R(v) of a node v ∈ V (S)
	// is empty if v has no left sibling or its left sibling w is
	// not stabbed. Otherwise, R(v) is the path from w to the rightmost
	// stabbed node in the subtree of w in S.
	var rightPath []*IntervalNode
	var nextNode *IntervalNode
	for nextNode = v.LeftSib; nextNode != nil; nextNode = nextNode.LeftSib {
		// If the left sibling stabs the query, then output, otherwise break
		if nextNode.GetChromStart() <= q.GetChromEnd() &&
			q.GetChromStart() <= nextNode.GetChromEnd() {
			rightPath = append(rightPath, nextNode)
		} else {
			break
		}
	}

	for _, node := range rightPath {
		traverseTree(node, stack, q, smaller)
	}
}

func stabbingQuery(q Interval, start *IntervalNode, smaller smallerIntervals) []Interval {
	stack := make([]Interval, 0)

	if start == nil {
		return stack
	}

	var pathToRoot []*IntervalNode
	var currNode *IntervalNode
	for currNode = start.Parent; currNode.Interval != nil; currNode = currNode.Parent {
		pathToRoot = append(pathToRoot, currNode)
	}

	for _, currNode = range pathToRoot {
		traverseTree(currNode, &stack, q, smaller)
	}

	return stack
}
