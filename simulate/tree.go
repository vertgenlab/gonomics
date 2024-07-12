package simulate

import (
	"fmt"
	"log"
	"math/rand"

	"github.com/vertgenlab/gonomics/expandedTree"
	"github.com/vertgenlab/gonomics/numbers"
)

// ETree produces a phylogenetic tree with a user-specified number of nodes and random gamma-distributed branch lengths,
// determined by a user-specific alpha and beta parameter.
func ETree(numNodes int, gammaAlpha float64, gammaBeta float64) *expandedTree.ETree {
	if numNodes%2 != 1 || numNodes < 0 {
		log.Fatalf("Error: Expecting a positive odd number of target nodes. Found: %v.\n", numNodes)
	}
	root := &expandedTree.ETree{
		Name: "root",
		Up:   nil,
	}
	var leaves = make([]*expandedTree.ETree, 0)
	generateChildNodes(root, gammaAlpha, gammaBeta, numNodes-1, leaves)
	return root
}

// generateChildNodes is a helper function of ETree which generates child nodes for a given node Up, recursively.
// At each recursive call, a random leaf is selected, and two child nodes are added to it. Branch lengths are
// randomly distributed based on a user-specified gamma distribution, which has two parameters alpha and beta.
// this process continues until the numNodesToAdd falls below 2.
func generateChildNodes(Up *expandedTree.ETree, gammaAlpha float64, gammaBeta float64, numNodesToAdd int, leaves []*expandedTree.ETree) {
	if numNodesToAdd < 2 {
		return
	}
	seed := rand.New(rand.NewSource(0))
	currBranchLength, _ := numbers.RandGamma(gammaAlpha, gammaBeta, seed)
	newLeftChild := &expandedTree.ETree{Name: fmt.Sprintf("Child_%v", numNodesToAdd), Up: Up, BranchLength: currBranchLength}
	Up.Left = newLeftChild
	currBranchLength, _ = numbers.RandGamma(gammaAlpha, gammaBeta, seed)
	newRightChild := &expandedTree.ETree{Name: fmt.Sprintf("Child_%v", numNodesToAdd-1), Up: Up, BranchLength: currBranchLength}
	Up.Right = newRightChild
	leaves = append(leaves, newLeftChild, newRightChild)
	nextUpIndex := rand.Intn(len(leaves))
	nextUp := leaves[nextUpIndex]
	leaves = append(leaves[:nextUpIndex], leaves[nextUpIndex+1:]...) //remove nextUp from the leaves slice
	generateChildNodes(nextUp, gammaAlpha, gammaBeta, numNodesToAdd-2, leaves)
}
