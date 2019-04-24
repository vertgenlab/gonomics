package reconstruct

import (
	"fmt"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/sophie/tree_newick"
)

//final function to run
func Reconstruct(root *tree_newick.NTree, filename_output string) {
	leaf := Get_leaf(root)
	branches := Get_branch(root)
	length := len(leaf[0].Fasta.Seq)

	for i := 0; i < length; i++ {
		//set up tree
		Set_state(root, i)
		Postorder(root)
		// loop nodes to get probable base at that site and node
		Loop_nodes(root)
	}
	var fastas []*fasta.Fasta

	for i := 0; i < len(branches); i++ {
		fastas = append(fastas, branches[i].Fasta)
	}

	var fas []fasta.Fasta
	//loop sites
	for i := 0; i < len(fastas); i++ {
		fas = append(fas, *fastas[i])
	}

	fasta.Write(filename_output, fas)

}

//test accuracy of the reconstruction compared to the simulation
func Accuracy(sim_filename string, rec_filename string) float64 {
	tot := 0.0
	sim, ers := fasta.Read(sim_filename)
	rec, err := fasta.Read(rec_filename)
	if ers != err {
	}
	des := "descendents_" + sim_filename
	sim_leafs, ersl := fasta.Read(des)
	if ersl != nil {
	}
	for i := 0; i < len(sim); i++ {
		for j := 0; j < len(sim_leafs); j++ {
			if sim[i].Name == sim_leafs[j].Name {
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

//get the nodes of the entire tree in a slice
func Get_tree(node *tree_newick.NTree) []*tree_newick.NTree {
	var branch []*tree_newick.NTree
	branch = append(branch, node)
	if node.Right != nil {
		b := Get_tree(node.Right)
		branch = append(branch, b...)
	}
	if node.Left != nil {
		a := Get_tree(node.Left)
		branch = append(branch, a...)
	}
	return branch
}

//get the interior nodes of the tree
func Get_branch(node *tree_newick.NTree) []*tree_newick.NTree {
	var branch []*tree_newick.NTree
	if node.Left != nil && node.Right != nil {
		branch = append(branch, node)
		a := Get_branch(node.Left)
		b := Get_branch(node.Right)
		branch = append(branch, a...)
		branch = append(branch, b...)
	}
	return branch
}

//get the leafs of the tree
func Get_leaf(node *tree_newick.NTree) []*tree_newick.NTree {
	var leaf []*tree_newick.NTree
	if node.Left != nil && node.Right != nil {
		a := Get_leaf(node.Left)
		b := Get_leaf(node.Right)
		leaf = append(leaf, a...)
		leaf = append(leaf, b...)
	}
	if node.Left == nil && node.Right == nil {
		leaf = append(leaf, node)
	}
	return leaf
}

//set the state of the tree given the Fasta and the position
func Set_state(node *tree_newick.NTree, pos int) {
	leaf := Get_leaf(node)
	var leaf_names []string
	for i := 0; i < len(leaf); i++ {
		leaf_names = append(leaf_names, leaf[i].Name)
	}
	node.State = 4
	node.Stored = []float64{0, 0, 0, 0}
	for i := 0; i < len(leaf_names); i++ {
		if node.Name == leaf_names[i] {
			if len(node.Fasta.Seq) > pos {
				node.State = int(node.Fasta.Seq[pos])
				node.Stored = []float64{0, 0, 0, 0}
				node.Stored[node.State] = 1
			}
		}
	}

	if node.Left != nil {
		Set_state(node.Left, pos)
	}
	if node.Right != nil {
		Set_state(node.Right, pos)
	}
}

//set up the tree with memory of the nodes below it using the set state
func Postorder(node *tree_newick.NTree) {
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
func Bubble_up(node *tree_newick.NTree, prev_node *tree_newick.NTree, scrap []float64) {
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
func Fix_fc(root *tree_newick.NTree, node *tree_newick.NTree) []float64 {
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
func Loop_nodes(root *tree_newick.NTree) {
	leafs := Get_leaf(root)
	branches := Get_branch(root)
	for j := 0; j < len(leafs[0].Fasta.Seq); j++ {
		for i := 0; i < len(leafs); i++ {
			leafs[i].State = int(leafs[i].Fasta.Seq[j])
		}
	}
	for j := 0; j < len(branches); j++ {
		//fix each interior node and get the most probable base
		fix := Fix_fc(root, branches[j])
		yhat := Yhat(fix)
		branches[j].Fasta.Seq = append(branches[j].Fasta.Seq, []dna.Base{dna.Base(yhat)}...)
	}
}
