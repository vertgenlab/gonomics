package simpleGraph

import (
	"sync"
	"fmt"
)

type Stack struct {
    Nodes []Node
    lock  sync.RWMutex
}

func NewStack(s *Stack) *Stack {
	s.Nodes = []Node{}
	return s
}

// Push adds an Item to the top of the stack
func Push(s *Stack, u Node) {
	s.lock.Lock()
	s.Nodes = append(s.Nodes, u)
	s.lock.Unlock()
}

// Pop removes an Item from the top of the stack
func Pop(s *Stack) *Node {
	s.lock.Lock()
	nodes := s.Nodes[len(s.Nodes)-1]
	s.Nodes = s.Nodes[0 : len(s.Nodes)-1]
	s.lock.Unlock()
	return &nodes
}

func (g *SimpleGraph) DFS(f func(*Node)) {
	g.lock.RLock()
	var st Stack
	var node *Node
	s := NewStack(&st)
	//add first node to the stack
	//n := g.Nodes[0]
	Push(s, *g.Nodes[0])
	visited := make(map[*Node]bool)

	var near []*Edge
	var i int
	var j *Node

	for {
		if s.IsEmpty() {
			break
		}
		node = Pop(s)
		visited[node] = true
		near = g.Edges[node]
		for i = 0; i < len(near); i++ {
			j = near[i].Next
			if !visited[j] {
				Push(s, *j)
				visited[j] = true
			}
		}
		if f != nil {
			f(node)
		}
	}
	g.lock.RUnlock()
}

func (s *Stack) IsEmpty() bool {
    s.lock.RLock()
    defer s.lock.RUnlock()
    return len(s.Nodes) == 0
}

func dfs(gg *SimpleGraph, node *Node, fn func (*Node)) {
    dfs_recur(gg, node, map[*Node]bool{}, fn)
}

func dfs_recur(gg *SimpleGraph, node *Node, v map[*Node]bool, fn func (*Node)) {

	 
    v[node] = true
    fn(node)
    for _, n := range gg.Edges[node] {
        if _, ok := v[n.Next]; !ok {
            dfs_recur(gg, n.Next, v, fn)
            fmt.Println(v)
        }
    }
}
