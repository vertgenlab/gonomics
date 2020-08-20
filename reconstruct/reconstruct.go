package reconstruct

import (
	"fmt"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/expandedTree"
	"github.com/vertgenlab/gonomics/fasta"
	"log"
)
//TODO: reconstruct never called
//final function to run
func Reconstruct(root *expandedTree.ETree) []*fasta.Fasta {
	leaves := expandedTree.GetLeaves(root)
	branches := expandedTree.GetBranch(root)
	var treeFastas []*fasta.Fasta

	seqLength := len(leaves[0].Fasta.Seq)
	for i := 0; i < seqLength; i++ {
		//set up tree
		SetLeafState(root, i)
		SetInternalState(root)
		// loop nodes to get probable base at that site and node
		LoopNodes(root)
	}
	for j := 0; j < len(branches); j++ {
		treeFastas = append(treeFastas, branches[j].Fasta)
	}

	return treeFastas

}

//test accuracy of the reconstruction compared to the simulation
func Accuracy(simFilename string, recFilename string) float64 {
	tot := 0.0 //why is this here when it could just be num?
	sim := fasta.Read(simFilename)
	rec := fasta.Read(recFilename)
	des := "descendents_" + simFilename
	simLeaves := fasta.Read(des)
	for i := 0; i < len(sim); i++ {
		for j := 0; j < len(simLeaves); j++ {
			if sim[i].Name == simLeaves[j].Name {
				sim = append(sim[:i], sim[i+1:]...)//this is adding the next record which hasn't been checked to the records?
				//sim is already all of the simulation fastas, why are we appending them?
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
		p = t / 3 //why branch length divided by 3?
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

//set the state of the tree given the Fasta and the position (zero-based)
func SetLeafState(node *expandedTree.ETree, pos int) {
	node.State = 4 // starts as N
	node.Stored = []float64{0, 0, 0, 0}

	if node.Right == nil && node.Left == nil { //if we are at a leaf
		if len(node.Fasta.Seq) <= pos {
			log.Fatal("position specified is out of range of sequence \n")
		} else if len(node.Fasta.Seq) > pos {
			node.State = int(node.Fasta.Seq[pos]) //State = base at given position
			for i := 0; i < 4; i++ { //check to make sure all other values of Stored are zero
				if i == node.State {
					node.Stored[i] = 1
				} else {
					node.Stored[i] = 0
				}
			}
		}
	}
}
//better description: set up Stored list for each node in the tree with probability of each base
//set up the tree with memory of the nodes below it using the set state
func SetInternalState(node *expandedTree.ETree) {
	if node.Left != nil && node.Right != nil {
		SetInternalState(node.Left)
		SetInternalState(node.Right)
		for i := 0; i < 4; i++ {
			sum := 0.0
			for j := 0; j < 4; j++ {
				for k := 0; k < 4; k++ {
					sum = sum + Prob(i, j, node.Left.BranchLength)*Prob(i, k, node.Right.BranchLength)*node.Left.Stored[j]*node.Right.Stored[k]
				}//probability that base i becomes either base j or k time the probability stored in the []float64 of that base appearing
			}
			node.Stored[i] = sum
		}
	} else if node.Left == nil && node.Right == nil{ //if at a leaf, use SetState to determine Stored Values
		for bases := 0; bases < len(node.Fasta.Seq); bases++ {
		SetLeafState(node, bases)
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
					if prevNode.Up.Left == prevNode {
						sum = sum + Prob(i, j, node.Left.BranchLength)*Prob(i, k, node.Right.BranchLength)*scrap[k]*node.Right.Stored[j]
					} else if prevNode.Up.Right == prevNode {
						sum = sum + Prob(i, j, node.Left.BranchLength)*Prob(i, k, node.Right.BranchLength)*scrap[j]*node.Left.Stored[k]
					} else {
						sum = sum + Prob(i, j, node.Left.BranchLength)*Prob(i, k, node.Right.BranchLength)*node.Left.Stored[k]*node.Right.Stored[j]
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
			//Bubble up the tree using the memory of the previous nodes given the fixed node
			BubbleUp(node.Up, node, scrap)
			ans[i] = root.Scrap //this is setting it to nothing?
		} else if node.Up == nil {
			ans[i] = root.Stored[i]
		}
	}

	return ans
}

//loop over the nodes of the tree to fix each node and append the most probable base to the Fasta
func LoopNodes(root *expandedTree.ETree) { //what is fixing the nodes doing?
	leaves := expandedTree.GetLeaves(root)
	branches := expandedTree.GetBranch(root)
	for j := 0; j < len(leaves[0].Fasta.Seq); j++ {
		for i := 0; i < len(leaves); i++ {
			leaves[i].State = int(leaves[i].Fasta.Seq[j])
		}
	}
	for j := 0; j < len(branches); j++ {
		//fix each interior node and get the most probable base
		fix := FixFc(root, branches[j])
		yHat := Yhat(fix)
		branches[j].Fasta.Seq = append(branches[j].Fasta.Seq, []dna.Base{dna.Base(yHat)}...)
	}
}
