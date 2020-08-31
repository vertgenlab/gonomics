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
				sim = append(sim[:i], sim[i+1:]...)
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
		p = t / 3
	}
	return p
}

//take in probability of all 4 bases return integer value of the most likely base
func Yhat(r []float64) int {
	var n float64
	n = 0
	var pos int
	for p, v := range r {
		if v > n {
			n = v
			pos = p
		}
	}
	return pos
}

func allZero(r []float64) bool {
	for _, v := range r {
		if v != 0 {
			return false
		}
	}
	return true
}

//set up Stored list for each node in the tree with probability of each base
func SetState(node *expandedTree.ETree, position int) {
	if node.Left != nil && node.Right != nil {
		SetState(node.Left, position)
		SetState(node.Right, position)
		if allZero(node.Left.Stored) && allZero(node.Right.Stored) {
			log.Fatal("no Stored values passed to internal node")
		}

		for i := 0; i < 4; i++ {
			sum := 0.0
			for j := 0; j < 4; j++ {
				for k := 0; k < 4; k++ {
					sum = sum + Prob(i, j, node.Left.BranchLength)*node.Left.Stored[j]*Prob(i, k, node.Right.BranchLength)*node.Right.Stored[k]
				}
			} //branch length (transition probability) times the probability of a base appearing given the context of the lower tree nodes
			node.Stored[i] = sum
		}
	} else if node.Left != nil {
		SetState(node.Left, position)
		if node.Left.Stored == nil {
			log.Fatal("no Stored values passed to internal node, left branch")
		}
		for i := 0; i < 4; i++ {
			sum := 0.0
			for j := 0; j < 4; j++ {
				sum = sum + Prob(i, j, node.Left.BranchLength)*node.Left.Stored[j]
			} //branch length (transition probability) times the probability of a base appearing given the context of the lower tree nodes
			node.Stored[i] = sum
		}
	} else if node.Right != nil {
		SetState(node.Right, position)
		if node.Right.Stored == nil {
			log.Fatal("no Stored values passed to internal node, right branch")
		}
		for i := 0; i < 4; i++ {
			sum := 0.0
			for k := 0; k < 4; k++ {
				sum = sum + Prob(i, k, node.Right.BranchLength)*node.Right.Stored[k]
			}
			node.Stored[i] = sum
		}
	} else if node.Right == nil && node.Left == nil {
		node.State = 4 // starts as N
		node.Stored = []float64{0, 0, 0, 0}

		if len(node.Fasta.Seq) <= position {
			log.Fatal("position specified is out of range of sequence \n")
		} else if len(node.Fasta.Seq) > position {
			node.State = int(node.Fasta.Seq[position])
			for i := 0; i < 4; i++ {
				if i == node.State {
					node.Stored[i] = 1
				} else {
					node.Stored[i] = 0
				}
			}
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
					if prevNode == node.Left { //scrap is equal to one position of prevNode.Stored (Left or Right)
						sum = sum + Prob(i, j, node.Left.BranchLength)*Prob(i, k, node.Right.BranchLength)*scrap[j]*node.Right.Stored[k]
					} else if prevNode == node.Right {
						sum = sum + Prob(i, j, node.Left.BranchLength)*Prob(i, k, node.Right.BranchLength)*scrap[k]*node.Left.Stored[j]
					}
				} else if prevNode.Up == nil {
					sum = sum + Prob(i, j, node.Left.BranchLength)*Prob(i, k, node.Right.BranchLength)*node.Left.Stored[j]*node.Right.Stored[k]
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
		scrap := []float64{0, 0, 0, 0} //checking one base at a time each time you call BubbleUp
		scrap[i] = node.Stored[i]
		if node.Up != nil {
			//Bubble up the tree using the memory of the previous node in relation to changing position taking in probabilities of bases
			//(node will be BubbleUp prevNode and node.Up will be the node being operated on)
			BubbleUp(node.Up, node, scrap) //node becomes PrevNode and scrap is set to one value of prevNode.Stored in BubbleUp
			ans[i] = root.Scrap            //root.Stored has previously assigned values (SetInternalState), you want to use whatever is returned by BubbleUp instead
		} else if node.Up == nil {
			ans[i] = root.Stored[i]
		}
	}

	return ans
}

//called by reconstructSeq.go on each base of the modern (leaf) seq. Loop over the nodes of the tree to return most probable base to the Fasta
func LoopNodes(root *expandedTree.ETree, position int) {
	internalNodes := expandedTree.GetBranch(root) //does this logic look like it will get leaves as well? We may want to use GetLeaves, which looks like it gets everything
	SetState(root, position)
	for k := 0; k < len(internalNodes); k++ {
		fix := FixFc(root, internalNodes[k])
		yHat := Yhat(fix)
		internalNodes[k].Fasta.Seq = append(internalNodes[k].Fasta.Seq, []dna.Base{dna.Base(yHat)}...)
	}
}
