package simpleGraph

import (
	"github.com/edotau/simpleio"
	"github.com/vertgenlab/gonomics/common"
	"github.com/vertgenlab/gonomics/dna"
	"log"
	"sync"
)

type MatrixAln struct {
	m     [][]int64
	trace [][]byte
}

func NewDnaPool() sync.Pool {
	return sync.Pool{
		New: func() interface{} {
			return make([]dna.Base, 0, 150)
		},
	}
}

func RightDynamicAln(alpha []dna.Base, beta []dna.Base, scores [][]int64, gapPen int64, m [][]int64, trace [][]byte) (int64, []simpleio.ByteCigar, int, int, int, int) {
	//check if size of alpha is larger than m
	var currMax int64
	var maxI int
	var maxJ int
	var i, j, routeIdx int
	//setting up the first rows and columns
	//seting up the rest of the matrix
	for i = 0; i < len(alpha)+1; i++ {
		for j = 0; j < len(beta)+1; j++ {
			if i == 0 && j == 0 {
				m[i][j] = 0
			} else if i == 0 {
				m[i][j] = m[i][j-1] + gapPen
				trace[i][j] = 'I'
			} else if j == 0 {
				m[i][j] = m[i-1][j] + gapPen
				trace[i][j] = 'D'
			} else {
				m[i][j], trace[i][j] = simpleio.ByteMatrixTrace(m[i-1][j-1]+scores[alpha[i-1]][beta[j-1]], m[i][j-1]+gapPen, m[i-1][j]+gapPen)
			}
			if m[i][j] > currMax {
				currMax = m[i][j]
				maxI = i
				maxJ = j
			}
		}
	}
	var route []simpleio.ByteCigar = make([]simpleio.ByteCigar, 0, 1)
	//traceback starts in top corner
	curr := simpleio.ByteCigar{}
	for i, j, routeIdx = maxI, maxJ, 0; i > 0 || j > 0; {
		//if route[routeIdx].RunLength == 0 {
		if len(route) == 0 {
			curr = simpleio.ByteCigar{RunLen: 1, Op: trace[i][j]}
			route = append(route, curr)
		} else if route[routeIdx].Op == trace[i][j] {
			route[routeIdx].RunLen += 1
		} else {
			curr = simpleio.ByteCigar{RunLen: 1, Op: trace[i][j]}
			route = append(route, curr)
			routeIdx++
		}
		switch trace[i][j] {
		case 'M':
			i, j = i-1, j-1
		case 'I':
			j -= 1
		case 'D':
			i -= 1
		default:
			log.Fatalf("Error: unexpected traceback with %c\n", trace[i][j])
		}
	}
	simpleio.ReverseBytesCigar(route)
	return m[maxI][maxJ], route, 0, maxI, 0, maxJ
}

func LeftDynamicAln(alpha []dna.Base, beta []dna.Base, scores [][]int64, gapPen int64, m [][]int64, trace [][]byte) (int64, []simpleio.ByteCigar, int, int, int, int) {
	//check if size of alpha is larger than m
	var i, j, routeIdx int

	for i = 0; i < len(alpha)+1; i++ {
		m[i][0] = 0
	}
	for j = 0; j < len(beta)+1; j++ {
		m[0][j] = 0
	}
	for i = 1; i < len(alpha)+1; i++ {
		for j = 1; j < len(beta)+1; j++ {
			m[i][j], trace[i][j] = simpleio.ByteMatrixTrace(m[i-1][j-1]+scores[alpha[i-1]][beta[j-1]], m[i][j-1]+gapPen, m[i-1][j]+gapPen)
			if m[i][j] < 0 {
				m[i][j] = 0
			}
		}
	}
	var minI, minJ = len(alpha), len(beta)
	var route []simpleio.ByteCigar = make([]simpleio.ByteCigar, 0, 1)
	curr := simpleio.ByteCigar{}
	//traceback starts in top corner
	for i, j, routeIdx = len(alpha), len(beta), 0; m[i][j] > 0; {
		//if route[routeIdx].RunLength == 0 {
		if len(route) == 0 {
			curr = simpleio.ByteCigar{RunLen: 1, Op: trace[i][j]}
			route = append(route, curr)
		} else if route[routeIdx].Op == trace[i][j] {
			route[routeIdx].RunLen += 1
		} else {
			curr = simpleio.ByteCigar{RunLen: 1, Op: trace[i][j]}
			route = append(route, curr)
			routeIdx++
		}
		switch trace[i][j] {
		case 'M':
			i, j = i-1, j-1
		case 'I':
			j -= 1
		case 'D':
			i -= 1
		default:
			log.Fatalf("Error: unexpected traceback")
		}
		minI = i
		minJ = j
	}
	//TODO: double check if this is tracing back in the correct directions
	simpleio.ReverseBytesCigar(route)
	return m[len(alpha)][len(beta)], route, minI, len(alpha), minJ, len(beta)
}

func getRightBases(n *Node, extenion int, tStartPos int, headSeq []dna.Base) []dna.Base {

	var targetLen int = common.Min(len(headSeq)+len(n.Seq)-tStartPos, extenion)
	var basesToTake int = targetLen - len(headSeq)
	var s []dna.Base = make([]dna.Base, targetLen)

	copy(s[0:len(headSeq)], headSeq)
	copy(s[len(headSeq):targetLen], n.Seq[tStartPos:tStartPos+basesToTake])
	return s
}

func RightAlignTraversal(n *Node, seq []dna.Base, start int, currentPath []uint32, ext int, read []dna.Base, m [][]int64, trace [][]byte, dnaPool *sync.Pool) ([]simpleio.ByteCigar, int64, int, int, []uint32) {

	s := dnaPool.Get().([]dna.Base)
	s = getRightBases(n, ext, start, seq)
	path := make([]uint32, 0, len(currentPath))
	copy(path, currentPath)
	var bestTargetEnd, bestQueryEnd, targetEnd, queryEnd int
	var bestScore, score int64
	var bestAlignment, alignment []simpleio.ByteCigar
	var bestPath []uint32

	if len(seq) >= ext {
		log.Fatalf("Error: the length of DNA sequence in previous nodes should not be enough to satisfy the desired extenion.\n")
	}

	//var availableBases int = len(seq) + len(n.Seq) - start
	//var targetLength int = common.Min(availableBases, ext)
	//var basesToTake int = targetLength - len(seq)

	//var s []dna.Base = make([]dna.Base, targetLength)
	//copy(s[0:len(seq)], seq)
	//copy(s[len(seq):targetLength], n.Seq[start:start+basesToTake])

	if len(s)-len(seq) >= ext || len(n.Next) == 0 {
		score, alignment, _, targetEnd, _, queryEnd = RightDynamicAln(s, read, HumanChimpTwoScoreMatrix, -600, m, trace)
		s = s[:0]
		dnaPool.Put(s)
		return alignment, score, targetEnd + start, queryEnd, path
	} else {
		bestScore = -1
		for _, i := range n.Next {
			alignment, score, targetEnd, queryEnd, path = RightAlignTraversal(i.Dest, s, 0, path, ext, read, m, trace, dnaPool)
			if score > bestScore {
				bestScore = score
				bestAlignment = alignment
				bestTargetEnd = targetEnd
				bestQueryEnd = queryEnd
				bestPath = path
			}
		}
	}
	return bestAlignment, bestScore, bestTargetEnd + start, bestQueryEnd, bestPath
}

//var cpuprofile = flag.String("cpuprofile", "", "write cpu profile to `file`")
//var memprofile = flag.String("memprofile", "", "write memory profile to `file`")

// Get Target sequence to align to
func getLeftTargetBases(n *Node, extenion int, tEndPos int, headSeq []dna.Base) []dna.Base {
	var targetLen int = common.Min(len(headSeq)+tEndPos, extenion)
	var basesToTake int = targetLen - len(headSeq)
	ans := make([]dna.Base, targetLen)
	copy(ans[0:basesToTake], n.Seq[tEndPos-basesToTake:tEndPos])
	copy(ans[basesToTake:targetLen], headSeq)
	return ans
}

func LeftAlignTraversal(n *Node, seq []dna.Base, refEnd int, currentPath []uint32, ext int, read []dna.Base, m [][]int64, trace [][]byte, dnaPool *sync.Pool) ([]simpleio.ByteCigar, int64, int, int, []uint32) {
	s := dnaPool.Get().([]dna.Base)
	s = getLeftTargetBases(n, ext, refEnd, seq)

	path := make([]uint32, 0, len(currentPath))
	copy(path, currentPath)
	AddPath(path, n.Id)
	var bestQueryStart, queryStart, refStart, bestRefStart int
	var bestScore, score int64
	var bestAlignment, alignment []simpleio.ByteCigar
	var bestPath []uint32

	if len(s)-len(seq) >= ext || len(n.Next) == 0 {
		score, alignment, refStart, _, queryStart, _ = LeftDynamicAln(s, read, HumanChimpTwoScoreMatrix, -600, m, trace)
		simpleio.ReversePath(path)
		s = s[:0]
		dnaPool.Put(s)
		return alignment, score, refEnd - len(s) - len(seq) + refStart, queryStart, path
	} else {
		bestScore = -1

		for _, i := range n.Prev {
			alignment, score, refStart, queryStart, path = LeftAlignTraversal(i.Dest, s, len(i.Dest.Seq), path, ext, read, m, trace, dnaPool)
			if score > bestScore {
				bestScore = score
				bestAlignment = alignment
				bestRefStart = refStart
				bestQueryStart = queryStart
				bestPath = path
			}
		}
	}

	return bestAlignment, bestScore, refEnd - len(s) - len(seq) + bestRefStart, bestQueryStart, bestPath
}

/*
func  Neighbors(n *Node) []Node {
	var nodes []Node
	for _, e := range n.edges {
		if ef == nil || ef(e) {
			if a := e.Tail(); a.ID() == n.ID() {
				nodes = append(nodes, e.Head())
			} else {
				nodes = append(nodes, a)
			}
		}
	}
	return nodes
}*/
