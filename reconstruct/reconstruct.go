// Package reconstruct has tools for reconstructing ancient genomes from the genomes of extant species using newick trees and fasta files.
package reconstruct

import (
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/expandedTree"
	"github.com/vertgenlab/gonomics/fasta"
	"log"
)

// WriteTreeToFasta writes assigned sequences at all nodes to a fasta file.
func WriteTreeToFasta(tree *expandedTree.ETree, outFile string) {
	var fastas []fasta.Fasta
	nodes := expandedTree.GetTree(tree)

	for i := 0; i < len(nodes); i++ {
		fastas = append(fastas, *nodes[i].Fasta)
	}
	fasta.Write(outFile, fastas)
}

// WriteLeavesToFasta writes assigned sequences at leaf nodes to a fasta file.
func WriteLeavesToFasta(tree *expandedTree.ETree, leafFile string) {
	var leafFastas []fasta.Fasta
	nodes := expandedTree.GetLeaves(tree)

	for i := 0; i < len(nodes); i++ {
		leafFastas = append(leafFastas, *nodes[i].Fasta)
	}
	fasta.Write(leafFile, leafFastas)
}

// mutationProbability calculate probability of switching from one base (a) to another (b) over branch length (t).
func mutationProbability(a int, b int, t float64) float64 {
	switch {
	case a > 3 || b > 3:
		return 0
	case a == b:
		return 1 - t
	default:
		return t / 3
	}
}

// LikelihoodsToBase returns the index of the most likely base for each position of a sequence. That position refers to a specific base.
func LikelihoodsToBase(likelihoods []float64, nonBiasBaseThreshold float64, biasBase dna.Base, highestProbThreshold float64) dna.Base {
	var highestProb, nonBiasBaseProb, total float64 = 0, 0, 0
	var answer dna.Base = biasBase
	for p, v := range likelihoods {
		total += v
		if dna.Base(p) != biasBase {
			nonBiasBaseProb += v
		}
		if v > highestProb {
			highestProb = v
			answer = dna.Base(p)
		}
	}
	if highestProb/total < highestProbThreshold {
		return dna.N
	}
	if nonBiasBaseProb/total < nonBiasBaseThreshold {
		return biasBase
	}
	return answer
}

// allZero checks if all values in a slice = 0 and returns a bool.
func allZero(r []float64) bool {
	for _, v := range r {
		if v != 0 {
			return false
		}
	}
	return true
}

// SetState calculates initial probabilities for all bases of the sequences for all nodes of the tree.
func SetState(node *expandedTree.ETree, position int) {
	var currNodeCounter, currLeftCounter, currRightCounter int
	var currSum float64

	if node.Left != nil && node.Right != nil {
		SetState(node.Left, position)
		SetState(node.Right, position)
		if allZero(node.Left.Stored) && allZero(node.Right.Stored) {
			log.Fatalf("Error: no Stored values passed to internal node.\n")
		}

		for currNodeCounter = range node.Stored {
			currSum = 0.0
			for currLeftCounter = range node.Left.Stored {
				for currRightCounter = range node.Right.Stored {
					currSum += mutationProbability(currNodeCounter, currLeftCounter, node.Left.BranchLength) * node.Left.Stored[currLeftCounter] * mutationProbability(currNodeCounter, currRightCounter, node.Right.BranchLength) * node.Right.Stored[currRightCounter]
				}
			} //branch length (transition probability) times the probability of a base appearing given the context of the lower tree nodes
			node.Stored[currNodeCounter] = currSum
		}
	} else if node.Left != nil {
		SetState(node.Left, position)
		if node.Left.Stored == nil {
			log.Fatalf("Error: no Stored values passed to internal node, left branch.\n")
		}
		for currNodeCounter = range node.Stored {
			currSum = 0.0
			for currLeftCounter = range node.Left.Stored {
				currSum += mutationProbability(currNodeCounter, currLeftCounter, node.Left.BranchLength) * node.Left.Stored[currLeftCounter]
			} //branch length (transition probability) times the probability of a base appearing given the context of the lower tree nodes
			node.Stored[currNodeCounter] = currSum
		}
	} else if node.Right != nil {
		SetState(node.Right, position)
		if node.Right.Stored == nil {
			log.Fatalf("Error: no Stored values passed to internal node, right branch.\n")
		}
		for currNodeCounter = range node.Stored {
			currSum = 0.0
			for currRightCounter = range node.Right.Stored {
				currSum += mutationProbability(currNodeCounter, currRightCounter, node.Right.BranchLength) * node.Right.Stored[currRightCounter]
			}
			node.Stored[currNodeCounter] = currSum
		}
	} else if node.Right == nil && node.Left == nil {
		node.State = 4 // starts as N
		node.Stored = []float64{0, 0, 0, 0}

		if len(node.Fasta.Seq) <= position {
			log.Fatal("position specified is out of range of sequence \n")
		} else {
			node.State = int(node.Fasta.Seq[position])
			for currNodeCounter = range node.Stored {
				if node.Fasta.Seq[position] == dna.N || node.Fasta.Seq[position] == dna.Gap {
					node.Stored[currNodeCounter] = 0.25
				} else if currNodeCounter == node.State {
					node.Stored[currNodeCounter] = 1
				} else {
					node.Stored[currNodeCounter] = 0
				}
			}
		}
	}
}

// bubbleUp calculates the final probabilities of all states in every position of the sequence at each internal node.
// using the stored values (initial probabilities from SetState) bubbleUp recursively calculates the
// probability on a child node based on both of the descendents of its ancestor. If a child node is a left child,
// bubbleUp uses the parent and right child's sequence information in Stored to compute a final probability
// of each of the base states at the left child, then passes those new probabilities up the tree to the root.
func bubbleUp(node *expandedTree.ETree, prevNode *expandedTree.ETree, scrap []float64) {
	var currNodeCounter, currLeftCounter, currRightCounter int
	var currSum, tot float64
	scrapNew := make([]float64, len(node.Stored))
	for currNodeCounter = range node.Stored {
		currSum = 0
		for currLeftCounter = range node.Left.Stored {
			for currRightCounter = range node.Right.Stored {
				if prevNode.Up != nil {
					if prevNode == node.Left { //scrap is equal to one position of prevNode.Stored (Left or Right)
						currSum += mutationProbability(currNodeCounter, currLeftCounter, node.Left.BranchLength) * mutationProbability(currNodeCounter, currRightCounter, node.Right.BranchLength) * scrap[currLeftCounter] * node.Right.Stored[currRightCounter]
					} else if prevNode == node.Right {
						currSum += mutationProbability(currNodeCounter, currLeftCounter, node.Left.BranchLength) * mutationProbability(currNodeCounter, currRightCounter, node.Right.BranchLength) * scrap[currRightCounter] * node.Left.Stored[currLeftCounter]
					}
				} else if prevNode.Up == nil {
					currSum += mutationProbability(currNodeCounter, currLeftCounter, node.Left.BranchLength) * mutationProbability(currNodeCounter, currRightCounter, node.Right.BranchLength) * node.Left.Stored[currLeftCounter] * node.Right.Stored[currRightCounter]
				}
			}
		}
		scrapNew[currNodeCounter] = currSum
	}
	if node.Up != nil {
		bubbleUp(node.Up, node, scrapNew)
	} else if node.Up == nil {
		tot = scrapNew[0] + scrapNew[1] + scrapNew[2] + scrapNew[3]
		node.Scrap = tot
	}
}

// FixFc passes a node to BubbleUp so that final base probabilities can be calculated from
// the initial values calculated in SetState.
func FixFc(root *expandedTree.ETree, node *expandedTree.ETree) []float64 {
	var currNodeCounter int
	ans := make([]float64, len(node.Stored))

	for currNodeCounter = range node.Stored {
		scrap := []float64{0, 0, 0, 0} //checking one base at a time each time you call BubbleUp
		scrap[currNodeCounter] = node.Stored[currNodeCounter]
		if node.Up != nil {
			//node will be bubbleUp prevNode and node.Up will be the node being operated on
			bubbleUp(node.Up, node, scrap)    //node becomes PrevNode and scrap is set to one value of prevNode.Stored in BubbleUp
			ans[currNodeCounter] = root.Scrap //root.Stored has previously assigned values (SetInternalState), you want to use whatever is returned by BubbleUp instead
		} else if node.Up == nil {
			ans[currNodeCounter] = root.Stored[currNodeCounter]
		}
	}

	return ans
}

func BaseExistsAtNodes(root *expandedTree.ETree, pos int) {
	descendentBaseExists(root, pos)
	baseExistsRecursive(root, pos)
}

// baseExistsRecursive is a helper function of BaseExistsAtNodes that determines, for a particular node at a particular
// position, whether a base exists (A, C, G, T, N)
// or not (node.Fasta.Seq[position] == dna.Gap). This information is stored as a bool in the field node.BasePresent.
// This function recursively searches descendent nodes and updates the BasePresent field.
func baseExistsRecursive(node *expandedTree.ETree, pos int) {
	var count int = 0
	if node.Left != nil && node.Right == nil {
		log.Fatalf("Error: tree is not a well-formed binary tree. For node: %v, left is not nil but right is nil.\n", node.Name)
	}
	if node.Right != nil && node.Left == nil {
		log.Fatalf("Error: tree is not a well-formed binary tree. For node: %v, right is not nil but left is nil.\n", node.Name)
	}
	if node.Up != nil && node.Up.BasePresent {
		count++
	}
	if node.Left != nil && node.Left.DescendentBasePresent {
		count++
	}
	if node.Right != nil && node.Right.DescendentBasePresent {
		count++
	}
	node.BasePresent = count >= 2
	if node.Left != nil {
		baseExistsRecursive(node.Left, pos)
	}
	if node.Right != nil {
		baseExistsRecursive(node.Right, pos)
	}
}

// descendentBaseExists is a helper function of BaseExistsAtNodes that determines whether, for a given node,
// any descendent nodes contain a base.
// This bool output is stored in node.DescendentBasePresent.
// This function recursively updates this field for all descendent nodes of the input node.
func descendentBaseExists(node *expandedTree.ETree, pos int) {
	if node.Left != nil && node.Right == nil {
		log.Fatalf("Error: tree is not a well-formed binary tree. For node: %v, left is not nil but right is nil.\n", node.Name)
	}
	if node.Right != nil && node.Left == nil {
		log.Fatalf("Error: tree is not a well-formed binary tree. For node: %v, right is not nil but left is nil.\n", node.Name)
	}
	if node.Left == nil && node.Right == nil { //if node is a leaf
		// a descendent base is present at a leaf if the sequence at this position is not a gap.
		node.DescendentBasePresent = node.Fasta.Seq[pos] != dna.Gap
	} else { //if we are not a leaf
		if node.Left != nil {
			descendentBaseExists(node.Left, pos)
		}
		if node.Right != nil {
			descendentBaseExists(node.Right, pos)
		}
		node.DescendentBasePresent = node.Left.DescendentBasePresent || node.Right.DescendentBasePresent
	}
}

// LoopNodes performs ancestral sequence reconstruction for all nodes of an input tree, specified by an input root node,
// at a user-specified alignment position.
// Options: the user may specify a 'biasLeafName'. When specified, the reconstruction of this sequence's immediate ancestor
// will be biased towards that descendent. The degree of this bias is controlled by the option 'nonBiasBaseThreshold'.
// The user may also specify a 'highestProbThreshold'. If the program is uncertain about ancestral reconstruction for a particular
// node, this option will allow LoopNodes to return an 'N' for that node instead.
func LoopNodes(root *expandedTree.ETree, position int, biasLeafName string, nonBiasBaseThreshold float64, highestProbThreshold float64) {
	var fix []float64
	var biasBase, answerBase dna.Base
	var biasParentName string
	var biasLeafNode *expandedTree.ETree

	if biasLeafName != "" {
		biasLeafNode = expandedTree.FindNodeName(root, biasLeafName)
		if biasLeafNode == nil {
			log.Fatalf("Didn't find %s in tree.\n", biasLeafName)
		}
		if biasLeafNode.Up == nil {
			log.Fatalf("Error: Bias reconstruction node was specified as the root node.")
		}
		biasParentName = biasLeafNode.Up.Name
	}

	internalNodes := expandedTree.GetBranch(root)
	SetState(root, position)
	BaseExistsAtNodes(root, position)
	for k := range internalNodes {
		fix = FixFc(root, internalNodes[k])

		if internalNodes[k].BasePresent {
			if biasParentName != "" && internalNodes[k].Name == biasParentName {
				biasBase = biasLeafNode.Fasta.Seq[position]
				answerBase = LikelihoodsToBase(fix, nonBiasBaseThreshold, biasBase, highestProbThreshold) //biased estimate
			} else {
				answerBase = LikelihoodsToBase(fix, 0, dna.N, highestProbThreshold) //unbiased estimate
			}
		} else {
			answerBase = dna.Gap
		}

		internalNodes[k].Fasta.Seq = append(internalNodes[k].Fasta.Seq, answerBase)
	}
}
