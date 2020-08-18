package reconstruct

import (
	"fmt"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/expandedTree"
	"github.com/vertgenlab/gonomics/fasta"
	"log"
)
//reconstruct never called
//final function to run
func Reconstruct(root *expandedTree.ETree) []*fasta.Fasta {
	leaves := expandedTree.GetLeaves(root)
	branches := expandedTree.GetBranch(root)
	var treeFastas []*fasta.Fasta

	seqLength := len(leaves[0].Fasta.Seq)
	for i := 0; i < seqLength; i++ {
		//set up tree
		SetState(root, i)
		Postorder(root)
		// loop nodes to get probable base at that site and node
		Loop_nodes(root)
	}
	for j := 0; j < len(branches); j++ {
		treeFastas = append(treeFastas, branches[j].Fasta)
	}

	return treeFastas

}

//test accuracy of the reconstruction compared to the simulation
func Accuracy(simFilename string, recFilename string) float64 {
	tot := 0.0
	sim := fasta.Read(simFilename)
	rec := fasta.Read(recFilename)
	des := "descendents_" + simFilename
	simLeafs := fasta.Read(des)
	for i := 0; i < len(sim); i++ {
		for j := 0; j < len(simLeafs); j++ {
			if sim[i].Name == simLeafs[j].Name {
				sim = append(sim[:i], sim[i+1:]...)
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
func SetState(node *expandedTree.ETree, pos int) {
	node.State = 4
	node.Stored = []float64{0, 0, 0, 0} //everything starts as N

	if node.Right == nil && node.Left == nil { //if we are at a leaf
		if len(node.Fasta.Seq) <= pos {
			log.Fatal("position specified is out of range of sequence \n")
		} else if len(node.Fasta.Seq) > pos {
			node.State = int(node.Fasta.Seq[pos]) //State = base at given position
		//the position in the list which corresponds to the base of State becomes 1 (state =0 (A) then 0 position of stored becomes 1)
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

//set up the tree with memory of the nodes below it using the set state
func Postorder(node *expandedTree.ETree) {
	if node.Left != nil && node.Right != nil {
		Postorder(node.Left)
		Postorder(node.Right)
		for i := 0; i < 4; i++ {
			sum := 0.0
			for j := 0; j < 4; j++ {
				for k := 0; k < 4; k++ {
					sum = sum + Prob(i, j, node.Left.BranchLength)*Prob(i, k, node.Right.BranchLength)*node.Left.Stored[j]*node.Right.Stored[k]
				}
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
func Bubble_up(node *expandedTree.ETree, prev_node *expandedTree.ETree, scrap []float64) {
	tot := 0.0
	scrap_new := []float64{0, 0, 0, 0}
	for i := 0; i < 4; i++ {
		sum := 0.0
		for j := 0; j < 4; j++ {
			for k := 0; k < 4; k++ {
				if prev_node.Up != nil {
					switch {
					case prev_node.Up.Left == prev_node:
						sum = sum + Prob(i, j, node.Left.BranchLength)*Prob(i, k, node.Right.BranchLength)*scrap[k]*node.Right.Stored[j]
					case prev_node.Up.Right == prev_node:
						sum = sum + Prob(i, j, node.Left.BranchLength)*Prob(i, k, node.Right.BranchLength)*scrap[j]*node.Left.Stored[k]
					default:
						sum = sum + Prob(i, j, node.Left.BranchLength)*Prob(i, k, node.Right.BranchLength)*node.Left.Stored[k]*node.Right.Stored[j]
					}
				} else if prev_node.Up == nil {
					sum = sum + Prob(i, j, node.Left.BranchLength)*Prob(i, k, node.Right.BranchLength)*node.Left.Stored[k]*node.Right.Stored[j]
				}
			}
		}
		scrap_new[i] = sum
	}
	if node.Up != nil {
		Bubble_up(node.Up, node, scrap_new)
	} else if node.Up == nil {
		tot = scrap_new[0] + scrap_new[1] + scrap_new[2] + scrap_new[3]
		node.Scrap = tot
	}
}

//fix each node and return the probabilities for each base at that site
func Fix_fc(root *expandedTree.ETree, node *expandedTree.ETree) []float64 {
	ans := []float64{0, 0, 0, 0}

	for i := 0; i < 4; i++ {
		scrap := []float64{0, 0, 0, 0}
		scrap[i] = node.Stored[i]
		if node.Up != nil {
			//Bubble up the tree using the memory of the previous nodes given the fixed node
			Bubble_up(node.Up, node, scrap)
			ans[i] = root.Scrap
		} else if node.Up == nil {
			ans[i] = root.Stored[i]
		}
	}

	return ans
}

//loop over the nodes of the tree to fix each node and append the most probable base to the Fasta
func Loop_nodes(root *expandedTree.ETree) {
	leaves := expandedTree.GetLeaves(root)
	branches := expandedTree.GetBranch(root)
	for j := 0; j < len(leaves[0].Fasta.Seq); j++ {
		for i := 0; i < len(leaves); i++ {
			leaves[i].State = int(leaves[i].Fasta.Seq[j])
		}
	}
	for j := 0; j < len(branches); j++ {
		//fix each interior node and get the most probable base
		fix := Fix_fc(root, branches[j])
		yhat := Yhat(fix)
		branches[j].Fasta.Seq = append(branches[j].Fasta.Seq, []dna.Base{dna.Base(yhat)}...)
	}
}
