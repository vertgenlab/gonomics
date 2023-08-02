package tree

import (
	"math"
)

func numLeaves(node *Tree) int {
	if node == nil {
		return 0
	} else {
		return 1 + numLeaves(node.Left) + numLeaves(node.Right)
	}
}

// Height returns the maximum cumulative branch length from the input tree node to a leaf node.
func Height(node *Tree) float64 {
	if node == nil {
		return 0
	} else {
		return node.BranchLength + math.Max(Height(node.Left), Height(node.Right))
	}
}
