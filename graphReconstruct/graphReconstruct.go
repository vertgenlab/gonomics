package graphReconstruct

import (
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/dnaTwoBit"

	//"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/expandedTree"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/simpleGraph"
	"log"
)

type graphColumn struct {
	AlignId    int
	AlignNodes map[string][]*simpleGraph.Node
}

//returns the percentage accuracy by base returned by reconstruct of each node and of all nodes combined (usage in reconstruct_test.go)
func ReconAccuracy(simFilename string, reconFilename string) map[string]float64 {
	var allNodes string
	allNodes = "all Nodes"
	var found bool = false
	var total float64
	total = 0.0
	var mistakes float64
	sim := fasta.Read(simFilename)
	recon := fasta.Read(reconFilename)

	answer := make(map[string]float64)

	for i := 0; i < len(sim); i++ {
		mistakes = 0.0
		found = false
		for j := 0; j < len(recon); j++ {
			if sim[i].Name == recon[j].Name {
				found = true
				//DEBUG: log.Printf("\n%s \n%s \n", dna.BasesToString(sim[i].Seq), dna.BasesToString(recon[j].Seq))
				for k := 0; k < len(sim[0].Seq); k++ {
					if sim[i].Seq[k] != recon[j].Seq[k] {
						mistakes = mistakes + 1
					}
				}
			}
		}
		if found == false {
			log.Fatal("Did not find all simulated sequences in reconstructed fasta.")
		}
		accuracy := mistakes / float64(len(sim[i].Seq)) * 100.0
		//DEBUG: fmt.Printf("tot: %f, len(sim): %f, len(sim[0].Seq): %f \n", tot, float64(len(sim)), float64(len(sim[0].Seq)))
		acc := 100 - accuracy
		answer[sim[i].Name] = acc
		total = total + mistakes
	}
	accuracy := total / (float64(len(sim)) * float64(len(sim[0].Seq))) * 100.0
	//DEBUG: fmt.Printf("tot: %f, len(sim): %f, len(sim[0].Seq): %f \n", tot, float64(len(sim)), float64(len(sim[0].Seq)))
	acc := 100 - accuracy
	answer[allNodes] = acc
	return answer
}

//write assigned sequences at all nodes to a fasta file
func WriteTreeToFasta(tree *expandedTree.ETree, outFile string) {
	var fastas []*fasta.Fasta
	nodes := expandedTree.GetTree(tree)

	for i := 0; i < len(nodes); i++ {
		fastas = append(fastas, nodes[i].Fasta)
	}
	fasta.Write(outFile, fastas)
}

//write assigned sequences at leaf nodes to a fasta file
func WriteLeavesToFasta(tree *expandedTree.ETree, leafFile string) {
	var leafFastas []*fasta.Fasta
	nodes := expandedTree.GetLeaves(tree)

	for i := 0; i < len(nodes); i++ {
		leafFastas = append(leafFastas, nodes[i].Fasta)
	}
	fasta.Write(leafFile, leafFastas)
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

func GraphProb(a float32, b float64) float64 {
	//calc probability that the edge.Prob is unchanged in parent
	answer := float64(a) * (1 - b)
	return answer
}

//set up Stored list for each node in the tree with probability of each base
func SetState(node *expandedTree.ETree, position int, graph *simpleGraph.SimpleGraph) *simpleGraph.SimpleGraph {
	for gn := 0; gn < len(graph.Nodes); gn++ {
		currNode := graph.Nodes[gn]
		incoming := currNode.Prev
		sum := 0.0
		if node.Left != nil && node.Right != nil {
			SetState(node.Left, position, graph)
			SetState(node.Right, position, graph)
			//outgoing := currNode.Next
			for i := 0; i < len(incoming); i++ {
				//only want a single edge calculation in sum
				sum = GraphProb(incoming[i].Prob, node.Left.BranchLength) * GraphProb(incoming[i].Prob, node.Right.BranchLength)
				//the probability that this node was passed through in node.r/l and the prob that edge.prob has remained the same (BranchLength)
				incoming[i].Prob = float32(sum)
				//updates existing edge (either this node's graph has been copied from below, and needs to have edges changed
				//or the graph is a single graph for the whole tree that gets passed up to the parent and the edges need to be updated)
			}
			//for i := 0; i < 4; i++ {
			//	sum := 0.0
			//	for j := 0; j < 4; j++ {
			//		for k := 0; k < 4; k++ {
			//			sum = sum + Prob(i, j, node.Left.BranchLength)*node.Left.Stored[j]*Prob(i, k, node.Right.BranchLength)*node.Right.Stored[k]
			//		}
			//	}
			//	node.Stored[i] = sum
			//}
		} else if node.Left != nil {
			SetState(node.Left, position, graph)
			//if node.Left.Stored == nil {
			//	log.Fatal("no Stored values passed to internal node, left branch")
			//}
			for i := 0; i < len(incoming); i++ {
				//only want a single edge calculation in sum
				sum = GraphProb(incoming[i].Prob, node.Left.BranchLength) * GraphProb(incoming[i].Prob, node.Right.BranchLength)
				//missing probs needs to be derived from the same EDGE
				//the probability that this node was passed through in node.r/l and the modifier that changes that prob
				incoming[i].Prob = float32(sum)
			}
		} else if node.Right != nil {
			SetState(node.Right, position, graph)
			//if node.Right.Stored == nil {
			//	log.Fatal("no Stored values passed to internal node, right branch")
			//}
			for i := 0; i < 4; i++ {
				sum := 0.0
				for k := 0; k < 4; k++ {
					sum = sum + Prob(i, k, node.Right.BranchLength)*node.Right.Stored[k]
				}
				node.Stored[i] = sum
			}
		} else if node.Right == nil && node.Left == nil {
			//node.State = 4 // starts as N
			//node.Stored = []float64{0, 0, 0, 0}
			if len(graph.Nodes) <= gn {
				log.Fatal("node is out of range")
			} //not sure that i need an else (i don't need to use set or stored yet or maybe at all. Don't need to store the base,
			// should maybe return a graph and then that can get fed to GraphRecon to create a likely path

			//if len(node.Fasta.Seq) <= position {
			//	log.Fatal("position specified is out of range of sequence \n")
			//} else if len(node.Fasta.Seq) > position {
			//	node.State = int(node.Fasta.Seq[position])
			//	for i := 0; i < 4; i++ {
			//		if i == node.State {
			//			node.Stored[i] = 1
			//		} else {
			//			node.Stored[i] = 0
			//		}
			//	}
			//}
		}
	}
	return graph
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
//func LoopNodes(root *expandedTree.ETree, position int) {
//	internalNodes := expandedTree.GetBranch(root)
//	SetState(root, position)
//	for k := 0; k < len(internalNodes); k++ {
//		fix := FixFc(root, internalNodes[k])
//		yHat := Yhat(fix)
//		internalNodes[k].Fasta.Seq = append(internalNodes[k].Fasta.Seq, []dna.Base{dna.Base(yHat)}...)
//	}
//}

//call for GraphRecon will be a lot like loopNodes call but instead of calling it on each position it will be called for each graph column
//most likely seq = > (branchLength x Prev) checking all Prev in each node of nodeInfo?
func GraphRecon(root *expandedTree.ETree, column graphColumn) {
	//var rightPresent = false
	//var leftPresent = false
	//nodeInfo the string is every Seq that exists in this alignColumn and the bool is a dummy, these will all become nodes for this species
	//making own graph for this species, making all possible nodes with all possible seqs (i.e. make a node for every node that exists in the daughters in this column)
	//if root.Right != nil && root.Left != nil {
	//	for c := 0; c < len(column.AlignNodes); c++ {
	//		for s := 0; s < len(column.AlignNodes[c]); s++ {
	//			if column.AlignNodes[c][s].Name == root.Right.Name { //node names must contain the species name
	//				rightPresent = true
	//				for name, info := range nodeInfo {
	//					if name == column.AlignNodes[c][s].Name {
	//						info = append(info, column.AlignNodes[c][s])
	//						nodeInfo[column.AlignNodes[c][s].Name] = info
	//					}
	//				}
	//			} else if column.AlignNodes[c][s].Name == root.Left.Name {
	//				leftPresent = true
	//				for name, info := range nodeInfo {
	//					if name == column.AlignNodes[c][s].Name {
	//						info = append(info, column.AlignNodes[c][s])
	//						nodeInfo[column.AlignNodes[c][s].Name] = info
	//					}
	//				}
	//			}
	//		}
	//	}
	//}
	//to assign next and prev i can look for a next/prev with a dest that does not exist in my nodeInfo map or in this align column, and make that the next/prev of the new node
	var nodeInfo = make(map[string]bool)
	for _, nodes := range column.AlignNodes {
		for n := range nodes {
			var newSeq = true
			var stringSeq string
			var check bool //DEBUG: to make sure we actually check the sequence, rather than fail to enter the loop
			for seq, _ := range nodeInfo {
				check = true
				if dna.BasesToString(nodes[n].Seq) == seq {
					newSeq = false
					stringSeq = dna.BasesToString(nodes[n].Seq)
				}
			}
			if newSeq && check { //TODO: remove check after debug
				nodeInfo[stringSeq] = true
			}
		}
	}

	//create unique num for each node in a graph?
	var num uint32
	for seq, _ := range nodeInfo {
		num += 1
		var newNode *simpleGraph.Node
		newNode = &simpleGraph.Node{Id: num, Name: root.Name, Seq: dna.StringToBases(seq), SeqTwoBit: dnaTwoBit.NewTwoBit(dna.StringToBases(seq)), Next: nil, Prev: nil, Info: simpleGraph.Annotation{}}
		for name, info := range column.AlignNodes { //add a new node to the existing nodes for this species
			if name == root.Name {
				info = append(info, newNode)
				column.AlignNodes[root.Name] = info
			} else {
				column.AlignNodes[root.Name] = []*simpleGraph.Node{newNode}
			}
		}
	}
}

//TODO: build a graph from a species' nodes in all columns
//maybe this function can add to it after each column is processed by graphRecon and then the new nodes being added can have their ID's changed
//according to the context of the whole graph

//seqOfPath takes in a graph and a path specified by the Node IDs and returns the seq of the path through the graph
func seqOfPath(g *simpleGraph.SimpleGraph, path []uint32) []dna.Base {
	var seq []dna.Base
	for n := 0; n < len(g.Nodes); n++ {
		for p := 0; p < len(path); p++ {
			if g.Nodes[n].Id == path[p] {
				seq = append(seq, g.Nodes[n].Seq...)
			}
		}
	}
	return seq
}

//PathFinder takes a graph and returns the most likely path through that graph after checking all possible paths from the first node to the last
func PathFinder(g *simpleGraph.SimpleGraph) ([]uint32, float32) {
	var finalPath []uint32
	var finalProb float32
	var tempPath = make([]uint32, 0)

	for n := 0; n < len(g.Nodes); n++ {
		if g.Nodes[n].Id == 0 {
			finalProb, finalPath = bestPath(g.Nodes[n], 1, tempPath)
		}
	}
	return finalPath, finalProb
}

//bestPath is the helper function for PathFinder, and recursively traverses the graph depth first to determine the most likely path from start to finish
func bestPath(node *simpleGraph.Node, prevProb float32, prevPath []uint32) (prob float32, path []uint32) {
	var tempProb float32 = 0
	var tempPath []uint32
	var finalProb float32
	var finalPath []uint32

	tempPath = append(tempPath, prevPath...)
	tempPath = append(tempPath, node.Id)
	if len(node.Next) == 0 {
		return prevProb, tempPath
	}
	for i, _ := range node.Next {
		tempProb = node.Next[i].Prob * prevProb
		currentProb, currentPath := bestPath(node.Next[i].Dest, tempProb, tempPath)
		if currentProb > finalProb {
			finalProb = currentProb
			finalPath = currentPath
		}
	}
	return finalProb, finalPath
}
