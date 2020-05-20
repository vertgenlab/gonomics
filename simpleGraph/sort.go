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

	// TODO: traverse forwards and backwards to get indegree table for each contiguous graph
	// Initialize starting state for inDegreeTable
	for i := 0; i < len(nodes); i++ {
		inDegreeTable[nodes[i]] = len(nodes[i].Prev)
	}

	// Find the first node with inDegree == 0
	for node, inDegree = range inDegreeTable {
		updatedNodes = nil
		if inDegree == 0 {
			answer = append(answer, node.Id)
			delete(inDegreeTable, node)
			updateTable(inDegreeTable, node, &updatedNodes)
			for k := 0; k < len(updatedNodes); k++ {
				answer = append(answer, updatedNodes[k].Id)
				delete(inDegreeTable, updatedNodes[k])
				updateTable(inDegreeTable, updatedNodes[k], &updatedNodes)
			}
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

func breakNonContiguousGraph(g *SimpleGraph) [][]*Node {
	answer := make([][]*Node, 0)
	var contiguousGraph []*Node
	inDegreeTable := make(map[*Node]int)
	visited := make([]bool, len(g.Nodes))
	var inDegree int
	var node *Node

	for i := 0; i < len(g.Nodes); i++ {
		inDegreeTable[g.Nodes[i]] = len(g.Nodes[i].Prev)
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

	for i = 0; i < len(startNode.Prev); i ++ {
		if !visited[startNode.Prev[i].Dest.Id] {
			*answer = append(*answer, startNode.Prev[i].Dest)
			visited[startNode.Prev[i].Dest.Id] = true
			traceGraph(startNode.Prev[i].Dest, visited, answer)
		}
	}
}