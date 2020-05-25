package simpleGraph

func GetSortOrder(g *SimpleGraph) []uint32 {
	return breadthFirstSearch(g.Nodes)
}

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

func updateTable(inDegreeTable map[*Node]int, node *Node, updatedNodes *[]*Node) {
	for i := 0; i < len(node.Next); i++ {
		inDegreeTable[node.Next[i].Dest]--
		if inDegreeTable[node.Next[i].Dest] == 0 {
			*updatedNodes = append(*updatedNodes, node.Next[i].Dest)
		}
	}
}

// TODO: possible to order nodes while breaking discontiguous graphs???
// TODO: presort graph node IDs and incorporate into simpleGraph??
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
