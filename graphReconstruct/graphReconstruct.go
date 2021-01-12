package graphReconstruct

import (
	//"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/expandedTree"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/simpleGraph"
	"log"
)

type graphColumn struct {
	AlignId    int
	AlignNodes [][]*simpleGraph.Node
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

//func GraphRecon(root *expandedTree.ETree, position int, graph simpleGraph.SimpleGraph) {
//	internalNodes := expandedTree.GetBranch(root)
//	//SetState(root, position, graph) //don't want a root graph
//	//loop align columns
//	for k := 0; k < len(internalNodes); k++ {
//		path := PathFinder(graph)
//		//make seq a dna.Base and make the va above path
//		//for i := 0; i < len(seq); i++ {
//			//internalNodes[k].Fasta.Seq = append(internalNodes[k].Fasta.Seq, seq[i])
//		}
//	}
//}

func PathFinder(g *simpleGraph.SimpleGraph) ([]uint32, float32) {
	var finalPath []uint32
	var finalProb float32

	for n := 0; n < len(g.Nodes); n++ {
		in := len(g.Nodes[n].Prev)
		if in == 0 {
			finalProb, finalPath = bestPath(*g.Nodes[n])
		}
	}

	return finalPath, finalProb
}

func bestPath(node simpleGraph.Node) (prob float32, path []uint32) {
	currentNode := node
	var currentProb float32 = 0
	currentPath := make([]uint32, 0)

	if len(currentNode.Next) == 0 {
		var lastNode []uint32
		lastNode = append(lastNode, currentNode.Id)
		return 1, lastNode
	}
	for i, _ := range currentNode.Next {
		tempProb, tempPath := bestPath(currentNode)
		if currentNode.Next[i].Prob*tempProb > currentProb {
			currentProb = currentNode.Next[i].Prob * tempProb
			currentPath = append(currentPath, currentNode.Id)
			currentPath = append(currentPath, tempPath...)
		}
	}
	return currentProb, currentPath
}

////PathFinder finds the best path from start to finish through the graph
//func PathFinder(g *simpleGraph.SimpleGraph) (finalPath []*simpleGraph.Node, prob float32) {
//	var path map[*simpleGraph.Node]int
//	var bestEdge int
//	var bestPath []*simpleGraph.Node
//	//var pathOptions = make(map[int][]float64) //map of paths options that contain a slice of the probs stored at the edges that it passes through
//	var nextNode *simpleGraph.Node
//	//var same = true
//	var visited = make(map[uint32]bool)
//	var lastNode = false
//	//var pQ []uint32
//	//var lastNode uint32
//	//lastNode = findLastNode(g)
//
//	for i := 0; i < len(g.Nodes); i++ {
//		in := len(g.Nodes[i].Prev)
//		out := len(g.Nodes[i].Next)
//		if in == 0 {
//			path[g.Nodes[i]] = 0
//			visited[g.Nodes[i].Id] = true
//		} else if in == 1 && out == 1 {
//			path[g.Nodes[i]] = 0
//			visited[g.Nodes[i].Id] = true
//		}
//		if out > 1 {
//			nextNode, bestEdge = pFHelper(g, i)
//			path[nextNode] = bestEdge
//			visited[g.Nodes[i].Id] = true
//			//start with edge that has highest prob
//			//traverse choosing the highest prob edges to the end.
//			//then start at the beginning and calculate cost, with only the non-highest edges (maybe just loop through and ignore highest prob edge?)
//			//(maybe if node has > 1 out && edge.Dest == !visited OR visited is true but there's more than one edge going into dest node of at least one of the out edges)
//			//and take it until it hits the end, if overall prob is lower than the previous traversal, ignore it
//			//trigger exit when all nodes have been visited, but that might happen before all possibilities are calculated
//			//
//			//
//			//and add nodeID to pQ
//			//follow to lastNode using the highest prob edges
//			//go to pQ, check if node.Next.Dest has been visited (if it has the edge that has the highest prob has already been used)
//			//if it has not been visited, i = the position of the []Nodes that this node is, take the most likely path from there to the lastNode
//			//While adding nodes to the pQ if they aren't already there
//			//at end check length of map with visited against number of nodes in graph, if ==, then set all visited bool to true and do conditional return
//			//need to add previous steps of the path into the map key that will hold each path, so as pQ has node added to it, create a new key referring to it,
//			//copy over the existing path that is being built and find a way to carry the variable that corresponds to that key of the map
//			//CornerCase: two edges have equal prob from a node
//			//PathFinder(g)
//			//} else if in > 1 { //finds a node where multiple paths converge
//			//	returnPath := findBestPath(g, i)
//			//	for l := 0; l < len(returnPath); l++ {
//			//		for p := 0; p < len(path); p++ {
//			//			if path[p].Id == returnPath[l].Id {
//			//				same = true
//			//			} else if path[p].Id != returnPath[l].Id {
//			//				same = false
//			//			}
//			//		}
//			//		if same == false {
//			//			path = append(path, returnPath[l])
//			//		}
//			//	}
//			//
//		} else if out == 0 {
//			lastNode = true
//			path[g.Nodes[i]] = -1
//		} else {
//			log.Fatal("unknown number of edges leading to this node")
//		}
//	}
//	pathProb := pathCalc(path)
//	for node, _ := range path {
//		bestPath = append(bestPath, node)
//	}
//
//	return bestPath, pathProb
//}
//
//func pathCalc(path map[*simpleGraph.Node]int) float32 {
//	var product float32
//	product = 1
//
//	for node, edge := range path {
//		product = product * node.Next[edge].Prob
//	}
//	return product
//}
//
//func pFHelper(g *simpleGraph.SimpleGraph, n int) (node *simpleGraph.Node, edge int) {
//	var probableNode *simpleGraph.Node
//	var nextNode *simpleGraph.Node
//	var bestEdge float32
//	var bE int
//	out := g.Nodes[n].Next
//
//	for e := 0; e < len(out); e++ {
//		if out[e].Prob > bestEdge {
//			bestEdge = out[e].Prob
//			nextNode = out[e].Dest //how to dereference Dest in this situation
//			bE = e
//		}
//	}
//	probableNode = nextNode
//	return probableNode, bE
//}
//
////findBestPath finds the most likely path when there is a divergence in the graph where multiple paths are possible, over the entire length of the divergence
//func findBestPath(graph *simpleGraph.SimpleGraph, n int) []*simpleGraph.Node {
//	//var possiblePaths map[int][]simpleGraph.Node
//	var simplePath []*simpleGraph.Node
//	var divergence = false
//
//	for k := 0; k < n; k++ {
//		out := len(graph.Nodes[k].Next)
//		if out == 1 {
//			divergence = false
//		} else if out < 1 {
//			divergence = true
//		}
//		if divergence {
//			//possiblePaths = make(map[int][]*simpleGraph.Node, out)
//
//		}
//	}
//
//	return simplePath
//}
//
////findLastNode goes through all the nodes of a graph and finds a node that has no Next Edges
//func findLastNode(g *simpleGraph.SimpleGraph) uint32 {
//	var lastNode uint32
//	for i := 0; i < len(g.Nodes); i++ {
//		out := g.Nodes[i].Next
//		if out == nil {
//			lastNode = g.Nodes[i].Id
//		}
//	}
//	return lastNode
//}
