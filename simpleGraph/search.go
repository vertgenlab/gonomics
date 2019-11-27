package simpleGraph

import (
	"fmt"
	"github.com/vertgenlab/gonomics/cigar"
	"github.com/vertgenlab/gonomics/common"
	"github.com/vertgenlab/gonomics/dna"
	"log"
	"sync"
)

func GraphTraversalFwd(g *SimpleGraph, n *Node, seq []dna.Base, path []uint32, start int, ext int) {
	s := make([]dna.Base, len(seq)+len(n.Seq)-start)
	path = AddPath(n.Id, path)
	copy(s[0:len(seq)], seq)
	copy(s[len(seq):len(seq)+len(n.Seq)-start], n.Seq[start:])
	if len(s) >= ext {
		//score, alignment, lowRef, _, lowQuery, highQuery = SmithWaterman(ref[common.StringToInt64(seedBeds[beds].Chrom)].Seq[seedBeds[beds].ChromStart:seedBeds[beds].ChromEnd], read.Seq, HumanChimpTwoScoreMatrix, -600, m, trace)
		fmt.Printf("Sequence: %s\n", dna.BasesToString(s[:ext]))
		//fmt.Printf("Path is: %s\n", path[0:len(path)-1])
	} else if len(n.Next) == 0 && len(s) < ext {
		fmt.Printf("Sequence: %s\n", dna.BasesToString(s))
		//fmt.Printf("Path is: %s\n", path[0:len(path)-1])
	} else {
		for _, i := range n.Next {
			//if path == "" {
			//	path = i.Dest.Name
			//} else {
			//	path =  path + ":" + i.Dest.Name
			//}
			GraphTraversalFwd(g, i.Dest, s, path, 0, ext)
		}
	}
}

func ReverseGraphTraversal(n *Node, seq []dna.Base, path []uint32, start int, ext int64) ([]uint32, []dna.Base) {
	s := make([]dna.Base, len(seq)+start)
	copy(s[0:start], n.Seq[:start])
	copy(s[start:start+len(seq)], seq)

	path = AddPath(n.Id, path)

	if int64(len(s)) >= ext {
		//fmt.Printf("Sequence: %s", dna.BasesToString(s[len(s)-ext:len(s)]))
		//fmt.Printf("Path is: %s\n", path)
		//GraphTraversalFwd(g, n, s[len(s)-ext:len(s)], path, start, ext)

	} else if len(n.Prev) == 0 && int64(len(s)) < ext {
		//fmt.Printf("Sequence: %s", dna.BasesToString(s))
		//fmt.Printf("Path is: %s\n", path)
		//GraphTraversalFwd(g, n, s, path, start, ext)
	} else {
		for _, i := range n.Prev {
			//fmt.Printf("Previous node: %s\n", dna.BasesToString(i.Next.Seq))
			path, s = ReverseGraphTraversal(i.Dest, s, path, len(i.Dest.Seq), ext)
		}
	}
	return path, s
}

func AddPath(newPath uint32, allPaths []uint32) []uint32 {
	if allPaths == nil {
		allPaths = append(allPaths, newPath)
	} else if allPaths[len(allPaths)-1] == newPath {
		return allPaths
	} else {
		allPaths = append(allPaths, newPath)
	}

	return allPaths
}

func reversePath(alpha []uint32) {
	for i, j := 0, len(alpha)-1; i < j; i, j = i+1, j-1 {
		alpha[i], alpha[j] = alpha[j], alpha[i]
	}
}

func PathToString(allPaths []uint32, gg *SimpleGraph) string {
	var s string = ""
	//fmt.Printf("length of paths %d\n", len(allPaths))
	if allPaths == nil {
		return s
	} else {
		s += gg.Nodes[allPaths[0]].Name
		if len(allPaths) > 1 {
			for i := 1; i < len(allPaths); i++ {
				s += ":" + gg.Nodes[allPaths[i]].Name
			}
		}

	}
	return s
}

func AlignTraversalFwd(n *Node, seq []dna.Base, start int, currentPath []uint32, ext int, read []dna.Base, m [][]int64, trace [][]rune) ([]*cigar.Cigar, int64, int, []uint32) {
	currentPath = append(currentPath, n.Id)
	var bestQueryEnd, queryEnd int
	var bestScore, score int64
	var bestAlignment, alignment []*cigar.Cigar
	var path, bestPath []uint32

	if len(seq) >= ext {
		log.Fatalf("Error: the length of DNA sequence in previous nodes should not be enough to satisfy the desired extenion.\n")
	}
	var availableBases int = len(seq) + len(n.Seq) - start
	var targetLength int = common.Min(availableBases, ext)
	var basesToTake int = targetLength - len(seq)
	var s []dna.Base = make([]dna.Base, targetLength)
	//log.Printf("len(seq)=%d, len(n.Seq)=%d, start=%d, targetLength=%d, basesToTake=%d\n", len(seq), len(n.Seq), start, targetLength, basesToTake)
	copy(s[0:len(seq)], seq)
	copy(s[len(seq):targetLength], n.Seq[start:start+basesToTake])

	if availableBases >= ext || len(n.Next) == 0 {
		score, alignment, _, _, _, queryEnd = RightLocal(s, read, HumanChimpTwoScoreMatrix, -600, m, trace)
		return alignment, score, queryEnd, currentPath
	} else {
		bestScore = -1
		for _, i := range n.Next {
			tmpPath := make([]uint32, len(currentPath)) //TODO: should not need to copy path N times, but N-1 times
			copy(tmpPath, currentPath)
			alignment, score, queryEnd, path = AlignTraversalFwd(i.Dest, s, 0, tmpPath, ext, read, m, trace)
			if score > bestScore {
				bestScore = score
				bestAlignment = alignment
				bestQueryEnd = queryEnd
				bestPath = path
			}
		}
		return bestAlignment, bestScore, bestQueryEnd, bestPath
	}
}

func AlignReverseGraphTraversal(n *Node, seq []dna.Base, refEnd int, currentPath []uint32, ext int, read []dna.Base, m [][]int64, trace [][]rune) ([]*cigar.Cigar, int64, int, int, []uint32) {
	currentPath = append([]uint32{n.Id}, currentPath...)
	var bestQueryStart, queryStart, refStart, bestRefStart int
	var bestScore, score int64
	var bestAlignment, alignment []*cigar.Cigar
	var path, bestPath []uint32

	var availableBases int = len(seq) + refEnd
	var targetLength int = common.Min(availableBases, ext)
	var basesToTake int = targetLength - len(seq)
	var s []dna.Base = make([]dna.Base, targetLength)
	copy(s[0:basesToTake], n.Seq[refEnd-basesToTake:refEnd])
	copy(s[basesToTake:targetLength], seq)

	//log.Printf("left(reverse) alignment: seq1=%s, seq2=%s\n", dna.BasesToString(s), dna.BasesToString(read))
	if availableBases >= ext || len(n.Next) == 0 {
		//log.Printf("at leaf, about to align, path is:%v\n", currentPath)
		score, alignment, refStart, _, queryStart, _ = LeftLocal(s, read, HumanChimpTwoScoreMatrix, -600, m, trace)
		return alignment, score, refEnd - basesToTake + refStart, queryStart, currentPath
	} else {
		bestScore = -1
		for _, i := range n.Prev {
			tmp := make([]uint32, len(currentPath))
			copy(tmp, currentPath)
			alignment, score, refStart, queryStart, path = AlignReverseGraphTraversal(i.Dest, s, len(i.Dest.Seq), currentPath, ext, read, m, trace)
			if score > bestScore {
				bestScore = score
				bestAlignment = alignment
				bestRefStart = refStart
				bestQueryStart = queryStart
				bestPath = path
			}
		}
		return bestAlignment, bestScore, bestRefStart, bestQueryStart, bestPath
	}
}

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

func (s *Stack) IsEmpty() bool {
	s.lock.RLock()
	defer s.lock.RUnlock()
	return len(s.Nodes) == 0
}

/*
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
}*/
