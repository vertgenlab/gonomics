package tree

import (
	"math"
)

/*
func branchHeight(node *Tree) float64 {
	if node == nil {
		return 0
	} else if node.OnlyTopology == true {
		return 1
	} else {
		return node.BranchLength
	}
}
*/

func numLeaves(node *Tree) int {
	if node == nil {
		return 0
	} else {
		return 1 + numLeaves(node.Left) + numLeaves(node.Right)
	}
}

func Height(node *Tree) float64 {
	if node == nil {
		return 0
	} else {
		return node.BranchLength + math.Max(Height(node.Left), Height(node.Right))
	}
}
