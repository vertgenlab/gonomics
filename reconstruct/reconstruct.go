package reconstruct

import (
	"fmt"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/expandedTree"
	"github.com/vertgenlab/gonomics/fasta"
	"log"
)

//replace reconstruct here with a cmd, loopNodes at the bottom now does everything this used to
//func Reconstruct(root *expandedTree.ETree) []*fasta.Fasta {
//	//leaves := expandedTree.GetLeaves(root)
//	branches := expandedTree.GetBranch(root)
//	var treeFastas []*fasta.Fasta
//
//	LoopNodes(root)
//
//	//seqLength := len(leaves[0].Fasta.Seq)
//	//for i := 0; i < seqLength; i++ { //loop through every base for SetLeafState
//	//	//set up tree
//	//	SetLeafState(root, i)
//	//	SetInternalState(root)
//	//	// loop nodes to get probable base at that site and node
//	//	LoopNodes(root)
//	//}
//	for j := 0; j < len(branches); j++ {
//		treeFastas = append(treeFastas, branches[j].Fasta)
//	}
//
//	return treeFastas
//}

//test accuracy of the reconstruction compared to the simulation
func Accuracy(simFilename string, recFilename string) float64 {
	tot := 0.0                     //why is this here when it could just be num?
	sim := fasta.Read(simFilename) //produced by simulate.printSeqForNodes, returned by simulate.Simulate
	rec := fasta.Read(recFilename)
	des := "descendents_" + simFilename //refers to file created by simulate.RemoveAncestors which has only leaf nodes labelled
	simLeaves := fasta.Read(des)
	for i := 0; i < len(sim); i++ {
		for j := 0; j < len(simLeaves); j++ {
			if sim[i].Name == simLeaves[j].Name { //if current fasta is a leaf fasta
				sim = append(sim[:i], sim[i+1:]...) //this is adding the next record?
				//sim is fastas of every node but the leaf nodes, so we are adding leaf fastas to the original simulation file
			}
		}
	}
	for i := 0; i < len(sim); i++ {
		num := 0.0

		for k := 0; k < len(sim[0].Seq); k++ {
			if sim[i].Seq[k] != rec[i].Seq[k] {
				num = num + 1
			}
		}
		tot = tot + num
	}
	accuracy := tot / (float64(len(sim)) * float64(len(sim[0].Seq))) * 100.0
	acc := 100 - accuracy
	fmt.Print("accuracy over all nodes= ", acc, "%", "\n")
	return acc
}

//calculate probability of switching from one base to another
func Prob(a int, b int, t float64) float64 {
	var p float64
	switch {
	case a > 3 || b > 3:
		p = 0
	case a == b:
		p = 1 - t
	default:
		p = t / 3 //chance there's a mutation divided among possible alternative outcomes
	}
	return p
}

//take in probability of all 4 bases return integer value of the most likely base
func Yhat(r []float64) int {
	var n float64
	var pos int
	for p, v := range r {
		if v > n {
			n = v
			pos = p
		}
	}
	return pos
}

//set the initial (leaf) state of the tree given the Fasta and the position (zero-based)
func SetLeafState(node *expandedTree.ETree, pos int) {
	node.State = 4 // starts as N
	node.Stored = []float64{0, 0, 0, 0}

	if node.Right == nil && node.Left == nil { //if we are at a leaf
		if len(node.Fasta.Seq) <= pos {
			log.Fatal("position specified is out of range of sequence \n")
		} else if len(node.Fasta.Seq) > pos {
			node.State = int(node.Fasta.Seq[pos]) //State = base at given position
			for i := 0; i < 4; i++ {              //check to make sure all other values of Stored are zero
				if i == node.State {
					node.Stored[i] = 1
				} else {
					node.Stored[i] = 0
				}
			}
		}
	} else if node.Right != nil {
		SetInternalState(node)
	} else if node.Left != nil {
		SetInternalState(node)
	}
}

//set up Stored list for each internal node in the tree with probability of each base
func SetInternalState(node *expandedTree.ETree) {
	if node.Left != nil {
		SetInternalState(node.Left)
		for i := 0; i < 4; i++ {
			sum := 0.0
			for j := 0; j < 4; j++ {
				sum = sum + Prob(i, j, node.Left.BranchLength)*node.Left.Stored[j]
			} //branch length (transition probability) times the probability of a base appearing given the context of the lower tree nodes
			node.Stored[i] = sum
		}
	} else if node.Right != nil {
		SetInternalState(node.Right)
		for i := 0; i < 4; i++ {
			sum := 0.0
			for k := 0; k < 4; k++ {
				sum = sum + Prob(i, k, node.Right.BranchLength)*node.Right.Stored[k]
			}
			node.Stored[i] = sum
		}
	} else {
		if node.State != 4 {
			node.Stored[node.State] = 1
		}
	}
}

//Bubble up the tree using the memory of the previous nodes
func BubbleUp(node *expandedTree.ETree, prevNode *expandedTree.ETree, scrap []float64) {
	tot := 0.0
	scrapNew := []float64{0, 0, 0, 0}
	for i := 0; i < 4; i++ {
		sum := 0.0
		for j := 0; j < 4; j++ {
			for k := 0; k < 4; k++ {
				if prevNode.Up != nil {
					if prevNode == node.Left {
						sum = sum + Prob(i, j, node.Left.BranchLength)*Prob(i, k, node.Right.BranchLength)*scrap[k]*node.Right.Stored[j]
					} else if prevNode == node.Right {
						sum = sum + Prob(i, j, node.Left.BranchLength)*Prob(i, k, node.Right.BranchLength)*scrap[j]*node.Left.Stored[k]
					}
				} else if prevNode.Up == nil {
					sum = sum + Prob(i, j, node.Left.BranchLength)*Prob(i, k, node.Right.BranchLength)*node.Left.Stored[k]*node.Right.Stored[j]
				}
			}
		}
		scrapNew[i] = sum
	}
	if node.Up != nil {
		BubbleUp(node.Up, node, scrapNew)
	} else if node.Up == nil {
		tot = scrapNew[0] + scrapNew[1] + scrapNew[2] + scrapNew[3]
		node.Scrap = tot
	}
}

//fix each node and return the probabilities for each base at that site
func FixFc(root *expandedTree.ETree, node *expandedTree.ETree) []float64 {
	ans := []float64{0, 0, 0, 0}

	for i := 0; i < 4; i++ {
		scrap := []float64{0, 0, 0, 0}
		scrap[i] = node.Stored[i] //not connecting this to root.Scrap?
		if node.Up != nil {
			//Bubble up the tree using the memory of the previous node in relation to changing position taking in probabilities of bases
			//(node will be BubbleUp prevNode and node.Up will be the node being operated on)
			BubbleUp(node.Up, node, scrap)
			ans[i] = root.Scrap
		} else if node.Up == nil {
			ans[i] = root.Stored[i]
		}
	}

	return ans
}

//loop over the nodes of the tree to fix each node and append the most probable base to the Fasta, return a slice of fastas for the whole tree
func LoopNodes(root *expandedTree.ETree) []*fasta.Fasta {
	leaves := expandedTree.GetLeaves(root)
	branches := expandedTree.GetBranch(root)
	var leafFastas []*fasta.Fasta
	var branchFastas []*fasta.Fasta
	var reconTree []*fasta.Fasta
	for i := 0; i < len(leaves[0].Fasta.Seq); i++ { //reconstruct a single position at every node in the tree and add it to that nodes Fasta
		for j := 0; j < len(leaves); j++ {
			SetLeafState(root, i)
			leaves[j].Fasta.Seq = append(leaves[j].Fasta.Seq, leaves[j].Fasta.Seq[i])
		}
		for k := 0; k < len(branches); k++ {
			SetInternalState(root)
			//fix each interior node and get the most probable base
			fix := FixFc(root, branches[k])
			yHat := Yhat(fix)
			branches[k].Fasta.Seq = append(branches[k].Fasta.Seq, []dna.Base{dna.Base(yHat)}...)
		}
	}
//now that each node has their reconstructed Fasta, create a slice of fastas for leaves and branches and return a single slice of fastas for the whole tree
	for l := 0; l < len(leaves); l++ {
		leafFastas = append(leafFastas, leaves[l].Fasta)
	}
	for lf := 0; lf < len(leafFastas); lf++ {
		reconTree = append(reconTree, leafFastas[lf])
	}
	for b := 0; b < len(branches); b++ {
		branchFastas = append(branchFastas, branches[b].Fasta)
	}
	for bf := 0; bf < len(branchFastas); bf++ {
		reconTree = append(reconTree, branchFastas[bf])
	}
	return reconTree
}
