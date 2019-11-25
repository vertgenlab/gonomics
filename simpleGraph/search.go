package simpleGraph

import (
	"fmt"
	"sync"
	//"log"
	"github.com/vertgenlab/gonomics/dna"
	//"github.com/vertgenlab/gonomics/fastq"
	"github.com/vertgenlab/gonomics/cigar"
	"github.com/vertgenlab/gonomics/sam"
)

func GraphTraversalFwd(g *SimpleGraph, n *Node, seq []dna.Base, path []int64, start int, ext int) {
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

func ReverseGraphTraversal(n *Node, seq []dna.Base, path []int64, start int, ext int64) ([]int64, []dna.Base) {
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

func AddPath(newPath int64, allPaths []int64) []int64 {
	if allPaths == nil {
		allPaths = append(allPaths, newPath)
	} else if allPaths[len(allPaths)-1] == newPath {
		return allPaths
	} else {
		allPaths = append(allPaths, newPath)
	}
	
	return allPaths
}

func reversePath(alpha []int64) {
	for i, j := 0, len(alpha)-1; i < j; i, j = i+1, j-1 {
		alpha[i], alpha[j] = alpha[j], alpha[i]
	}
}

func PathToString(allPaths []int64, gg *SimpleGraph) string {
	var s string = ""
	//fmt.Printf("length of paths %d\n", len(allPaths))
	if allPaths == nil {
		return s
	} else {
		s += gg.Nodes[allPaths[0]].Name
		if len(allPaths) > 1 {
			for i := 1; i < len(allPaths);i++ {
				s+= ":" + gg.Nodes[allPaths[i]].Name
			}
		}
		
	}
	return s
}

func AlignTraversalFwd(n *Node, seq []dna.Base, start int, bestPath []int64, ext int64, read []dna.Base, m [][]int64, trace [][]rune, bestCigar []*cigar.Cigar, bestScore int64, queryEnd int64) ([]*cigar.Cigar, int64, int64, []int64) {
	s := make([]dna.Base, len(seq)+len(n.Seq)-start)
	bestPath = AddPath(n.Id, bestPath)
	var score, maxJ int64

	var alignment []*cigar.Cigar
	copy(s[0:len(seq)], seq)
	copy(s[len(seq):len(seq)+len(n.Seq)-start], n.Seq[start:])

	if int64(len(s)) >= ext {
		score, alignment, _, _, _, maxJ = RightLocal(s[:ext], read, HumanChimpTwoScoreMatrix, -600, m, trace)
		if score > bestScore {
			bestScore = score
			bestCigar = alignment
			queryEnd = int64(len(read)) + maxJ
		}
		return bestCigar, bestScore, queryEnd, bestPath

	} else if len(n.Next) == 0 && int64(len(s)) < ext {

		score, alignment, _, _, _, maxJ = RightLocal(s, read, HumanChimpTwoScoreMatrix, -600, m, trace)
		if score > bestScore {
			bestScore = score
			bestCigar = alignment
			queryEnd = int64(len(read)) + maxJ
		}
		return bestCigar, bestScore, queryEnd, bestPath
	} else {
		for _, i := range n.Next {
			//bestPath += i.Dest.Name + ":"
			
			AlignTraversalFwd(i.Dest, s, 0, bestPath, ext, read, m, trace, bestCigar, bestScore, queryEnd)
		}
	}
	return bestCigar, bestScore, queryEnd, bestPath
}

func AlignReverseGraphTraversal(n *Node, seq []dna.Base, start uint64, bestPath []int64, ext int64, read []dna.Base, m [][]int64, trace [][]rune, currBest *sam.SamAln, bestScore int64, queryStart int64) (*sam.SamAln, int64, int64, []int64) {
	s := make([]dna.Base, uint64(len(seq))+start)
	//bestPath = AddPath(n.Name, bestPath)
	//if bestPath == "" {
	//	bestPath += n.Name
	//} else {
	//	bestPath = n.Name + ":" + bestPath
	//}
	bestPath = AddPath(n.Id, bestPath)
	var score, minJ int64 = 0, 0
	//var minJ int64
	//var lowRef, lowQuery, highQuery int64
	var alignment []*cigar.Cigar

	copy(s[0:start], n.Seq[:start])
	copy(s[start:start+uint64(len(seq))], seq)
	var refStart int64
	if int64(len(s)) >= ext {
		//fmt.Printf("Sequence: %s\n", dna.BasesToString(s[int64(len(s))-ext:len(s)]))
		score, alignment, refStart, _, minJ, _ = LeftLocal(s[int64(len(s))-ext:len(s)], read, HumanChimpTwoScoreMatrix, -600, m, trace)

		//fmt.Printf("Sequence: %s\n", dna.BasesToString(s[:ext]))
		if score > bestScore {
			reversePath(bestPath)
			bestScore = score
			currBest.Flag = 0
			//currBest.RName = bestPath[0]
			//log.Printf("Length of s: %d, ext: %d, refStart: %d", len(s), ext, refStart)
			currBest.Pos = int64(len(s)) - ext + refStart
			currBest.Cigar = alignment
			queryStart = int64(len(s)) - minJ
			//currBest.Extra = PathToString(path)
		}
		return currBest, bestScore, queryStart, bestPath

	} else if len(n.Prev) == 0 && int64(len(s)) < ext {
		//fmt.Printf("Sequence: %s\n", dna.BasesToString(s))

		score, alignment, refStart, _, minJ, _ = LeftLocal(s, read, HumanChimpTwoScoreMatrix, -600, m, trace)
		if score > bestScore {
			reversePath(bestPath)
			bestScore = score
			currBest.Flag = 0
			//currBest.RName = bestPath[0]
			currBest.Pos = int64(len(s)) - ext + refStart
			currBest.Cigar = alignment
			queryStart = int64(len(s)) - minJ
		}
		return currBest, bestScore, queryStart, bestPath
	} else {
		for _, i := range n.Prev {
			//bestPath = i.Dest.Name + ":" + bestPath
			//fmt.Printf("Previous node: %s\n", dna.BasesToString(i.Next.Seq))
			currBest, bestScore, queryStart, bestPath = AlignReverseGraphTraversal(i.Dest, s, uint64(len(i.Dest.Seq)), bestPath, ext, read, m, trace, currBest, bestScore, queryStart)
		}
	}
	return currBest, bestScore, queryStart, bestPath
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
