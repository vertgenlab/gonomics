package simpleGraph

import (
	"sync"
	"fmt"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/fastq"
	"github.com/vertgenlab/gonomics/sam"
	"github.com/vertgenlab/gonomics/cigar"
)

func GraphTraversalFwd(g *SimpleGraph, n *Node, seq []dna.Base, path string, start int, ext int) {
	s := make([]dna.Base, len(seq) + len(n.Seq)-start)
	path = path + n.Name + ":"
	copy(s[0:len(seq)], seq)
	copy(s[len(seq):len(seq)+len(n.Seq)-start], n.Seq[start:])
	if len(s) >= ext {
		//score, alignment, lowRef, _, lowQuery, highQuery = SmithWaterman(ref[common.StringToInt64(seedBeds[beds].Chrom)].Seq[seedBeds[beds].ChromStart:seedBeds[beds].ChromEnd], read.Seq, HumanChimpTwoScoreMatrix, -600, m, trace)
		fmt.Printf("Sequence: %s\n", dna.BasesToString(s[:ext]))
		fmt.Printf("Path is: %s\n", path[0:len(path)-1])
	} else if len(n.Next) == 0 && len(s) < ext {
		fmt.Printf("Sequence: %s\n", dna.BasesToString(s))
		fmt.Printf("Path is: %s\n", path[0:len(path)-1])
	} else {
		for _, i := range n.Next {
			
			GraphTraversalFwd(g, i.Next, s, path, 0, ext)
		}
	}
}

func ReverseGraphTraversal(g *SimpleGraph, n *Node, seq []dna.Base, path string, start int, ext int) (string, []dna.Base) {
	s := make([]dna.Base, len(seq) + start)
	copy(s[0:start], n.Seq[:start])


	copy(s[start:start+len(seq)], seq)
	if len(s) >= ext {
		//fmt.Printf("Sequence: %s", dna.BasesToString(s[len(s)-ext:len(s)]))
		//fmt.Printf("Path is: %s\n", path)
		//GraphTraversalFwd(g, n, s[len(s)-ext:len(s)], path, start, ext)
		
	} else if len(n.Prev) == 0 && len(s) < ext {
		//fmt.Printf("Sequence: %s", dna.BasesToString(s))
		//fmt.Printf("Path is: %s\n", path)
		//GraphTraversalFwd(g, n, s, path, start, ext)
	} else {
		for _,i := range n.Prev {
			//fmt.Printf("Previous node: %s\n", dna.BasesToString(i.Next.Seq))
			path = i.Next.Name + ":" + path
			path, s = ReverseGraphTraversal(g, i.Next, s, path, len(i.Next.Seq), ext)
		}
	}
	return path, s
}

func AlignTraversalFwd(g *SimpleGraph, n *Node, seq []dna.Base, start int, bestPath string, ext int, read fastq.Fastq, m [][]int64, trace [][]rune, currBest *sam.SamAln, bestScore int64) (*sam.SamAln, int64) {
	s := make([]dna.Base, len(seq) + len(n.Seq)-start)
	var score int64 = 0
	var lowRef int64
	//var lowRef, lowQuery, highQuery int64
	var alignment []*cigar.Cigar
	copy(s[0:len(seq)], seq)
	copy(s[len(seq):len(seq)+len(n.Seq)-start], n.Seq[start:])

	bestPath += n.Name + ":"
	
	if len(s) >= ext {
		score, alignment, lowRef, _, _, _ = SmithWaterman(s[:ext], read.Seq, HumanChimpTwoScoreMatrix, -600, m, trace)
		
		//fmt.Printf("Sequence: %s\n", dna.BasesToString(s[:ext]))
		if score > bestScore {
			bestScore = score
			currBest.Flag = 0
			currBest.RName = bestPath[0:len(bestPath)-1]
			//+1 to convert into one base sam record
			currBest.Pos = lowRef
			currBest.Cigar = alignment
			currBest.Seq = read.Seq
			currBest.Qual = string(read.Qual)
		}
	} else if len(n.Next) == 0 && len(s) < ext {
		//fmt.Printf("Sequence: %s\n", dna.BasesToString(s))
		score, alignment, lowRef, _, _, _ = SmithWaterman(s, read.Seq, HumanChimpTwoScoreMatrix, -600, m, trace)
		if score > bestScore {
			bestScore = score
			currBest.Flag = 0
			currBest.RName = bestPath[0:len(bestPath)-1]
			currBest.Pos = lowRef
			currBest.Cigar = alignment
			currBest.Seq = read.Seq
			currBest.Qual = string(read.Qual)
		}
	} else {
		for _, i := range n.Next {
			currBest, bestScore = AlignTraversalFwd(g, i.Next, s, 0, bestPath, ext, read, m, trace, currBest, bestScore)
		}
	}
	return currBest, bestScore
}

func AlignReverseGraphTraversal(g *SimpleGraph, n *Node, seq []dna.Base, start int, bestPath string, ext int, read fastq.Fastq, m [][]int64, trace [][]rune, currBest *sam.SamAln, bestScore int64) (*sam.SamAln, int64) {
	s := make([]dna.Base, len(seq) + start)
	bestPath += n.Name + ":"
	var score int64 = 0
	var lowRef int64
	//var lowRef, lowQuery, highQuery int64
	var alignment []*cigar.Cigar

	copy(s[0:start], n.Seq[:start])
	copy(s[start:start+len(seq)], seq)
	if len(s) >= ext {
		fmt.Printf("Sequence: %s\n", dna.BasesToString(s[len(s)-ext:len(s)]))
		score, alignment, lowRef, _, _, _ = SmithWaterman(s[len(s)-ext:len(s)], read.Seq, HumanChimpTwoScoreMatrix, -600, m, trace)
		
		//fmt.Printf("Sequence: %s\n", dna.BasesToString(s[:ext]))
		if score > bestScore {
			bestScore = score
			currBest.Flag = 0
			currBest.RName = bestPath[0:len(bestPath)-1]
			currBest.Pos = currBest.Pos+lowRef
			currBest.Cigar = alignment
			currBest.Seq = read.Seq
			currBest.Qual = string(read.Qual)
		}
	} else if len(n.Prev) == 0 && len(s) < ext {
		fmt.Printf("Sequence: %s\n", dna.BasesToString(s))

		score, alignment, lowRef, _, _, _ = SmithWaterman(s, read.Seq, HumanChimpTwoScoreMatrix, -600, m, trace)
		if score > bestScore {
			bestScore = score
			currBest.Flag = 0
			currBest.RName = bestPath[0:len(bestPath)-1]
			currBest.Pos = currBest.Pos+lowRef
			currBest.Cigar = alignment
			currBest.Seq = read.Seq
			currBest.Qual = string(read.Qual)
		}
	} else {
		for _,i := range n.Prev {
			//fmt.Printf("Previous node: %s\n", dna.BasesToString(i.Next.Seq))
			currBest, bestScore = AlignReverseGraphTraversal(g, i.Next, s, len(i.Next.Seq), bestPath, ext, read, m, trace, currBest, bestScore)
		}
	}
	return currBest, bestScore
}
/*
func dfs(gg *SimpleGraph, node *Node, fn func (*Node)) {
    dfsRecursion(gg, node, map[*Node]bool{}, fn)
}

func dfsRecursion(gg *SimpleGraph, node *Node, v map[*Node]bool, fn func (*Node)) {
    v[node] = true
    fn(node)
    for _, n := range gg.Edges[node] {
        if _, ok := v[n.Next]; !ok {
            dfsRecursion(gg, n.Next, v, fn)
        }
    }
}*/

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
