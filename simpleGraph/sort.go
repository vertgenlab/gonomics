package simpleGraph

// SortGraph will reorder nodes in a graph such that the order and Ids of the output graph are topologically sorted
func SortGraph(g *SimpleGraph) *SimpleGraph {
	answer := new(SimpleGraph)
	answer.Nodes = make([]*Node, len(g.Nodes))
	order := GetSortOrder(g)
	for sortedIdx, originalIdx := range order {
		answer.Nodes[sortedIdx] = g.Nodes[originalIdx]
		answer.Nodes[sortedIdx].Id = uint32(sortedIdx)

	}
	return answer
}

// GetSortOrder will perform a breadth first search (BFS) on a graph and return an output slice where output[sortedIdx] = originalIdx
func GetSortOrder(g *SimpleGraph) []uint32 {
	return breadthFirstSearch(g.Nodes)
}

// TODO: design function to get start positions only
// breadthFirstSearch performs a breadth first search on a graph and returns a slice correlating the sort order to the original order
func breadthFirstSearch(nodes []*Node) []uint32 {
	answer := make([]uint32, 0)
	var inDegree int
	var node *Node
	inDegreeTable := make(map[*Node]int)

	// Updated nodes is going to keep track of each node
	// which has had a change to it's inDegree
	// We will use this to loop through the graph
	// visiting each group of connected nodes in order
	// and by searching within each group with a
	// breadth-first approach
	updatedNodes := make([]*Node, 0)

	subGraphs := BreakNonContiguousGraph(nodes)

	// loop through each contiguous subGraph
	for _, nodeSet := range subGraphs {
		updatedNodes = nil
		inDegreeTable = make(map[*Node]int)
		for i := 0; i < len(nodeSet); i++ {
			inDegreeTable[nodeSet[i]] = len(nodeSet[i].Prev)
		}

		// Find all nodes that start with inDegree zero and add to updatedNodes
		for node, inDegree = range inDegreeTable {
			if inDegree == 0 {
				updatedNodes = append(updatedNodes, node)
			}
		}

		for k := 0; k < len(updatedNodes); k++ {
			answer = append(answer, updatedNodes[k].Id)
			delete(inDegreeTable, updatedNodes[k])
			updateTable(inDegreeTable, updatedNodes[k], &updatedNodes)
		}
	}
	return answer
}

// updateTable updates the table of node in degrees
func updateTable(inDegreeTable map[*Node]int, node *Node, updatedNodes *[]*Node) {
	for i := 0; i < len(node.Next); i++ {
		inDegreeTable[node.Next[i].Dest]--
		if inDegreeTable[node.Next[i].Dest] == 0 {
			*updatedNodes = append(*updatedNodes, node.Next[i].Dest)
		}
	}
}

// TODO: possible to order nodes while breaking discontiguous graphs???
// BreakNonContiguousGraph will return a slice of graphs ([]*Node) such that each graph in the slice is contiguous
func BreakNonContiguousGraph(g []*Node) [][]*Node {
	answer := make([][]*Node, 0)
	var contiguousGraph []*Node
	inDegreeTable := make(map[*Node]int)
	visited := make([]bool, len(g))
	var inDegree int
	var node *Node

	for i := 0; i < len(g); i++ {
		inDegreeTable[g[i]] = len(g[i].Prev)
	}

	for node, inDegree = range inDegreeTable {
		if inDegree == 0 && !visited[node.Id] {
			contiguousGraph = make([]*Node, 1)
			contiguousGraph[0] = node
			visited[node.Id] = true
			traceGraph(node, visited, &contiguousGraph)
			answer = append(answer, contiguousGraph)
		}
	}

	return answer
}

// traceGraph is a helper function that traverses graph and keeps track of which nodes have been visited
func traceGraph(startNode *Node, visited []bool, answer *[]*Node) {
	var i int = 0

	for i = 0; i < len(startNode.Next); i++ {
		if !visited[startNode.Next[i].Dest.Id] {
			*answer = append(*answer, startNode.Next[i].Dest)
			visited[startNode.Next[i].Dest.Id] = true
			traceGraph(startNode.Next[i].Dest, visited, answer)
		}
	}

	for i = 0; i < len(startNode.Prev); i++ {
		if !visited[startNode.Prev[i].Dest.Id] {
			*answer = append(*answer, startNode.Prev[i].Dest)
			visited[startNode.Prev[i].Dest.Id] = true
			traceGraph(startNode.Prev[i].Dest, visited, answer)
		}
	}
}
