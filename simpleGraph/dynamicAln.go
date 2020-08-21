package simpleGraph


import(
	"github.com/vertgenlab/gonomics/dna"
	"log"
)

func RightDynamicAln(alpha []dna.Base, beta []dna.Base, scores [][]int64, gapPen int64, m [][]int64, trace [][]CigarOp) (int64, []ByteCigar, int, int, int, int) {
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
				m[i][j], trace[i][j] = ByteMatrixTrace(m[i-1][j-1]+scores[alpha[i-1]][beta[j-1]], m[i][j-1]+gapPen, m[i-1][j]+gapPen)
			}
			if m[i][j] > currMax {
				currMax = m[i][j]
				maxI = i
				maxJ = j
			}
		}
	}
	var route []ByteCigar = make([]ByteCigar, 0, 1)
	//traceback starts in top corner
	curr := ByteCigar{}
	for i, j, routeIdx = maxI, maxJ, 0; i > 0 || j > 0; {
		//if route[routeIdx].RunLength == 0 {
		if len(route) == 0 {
			curr = ByteCigar{RunLen: 1, Op: trace[i][j]}
			route = append(route, curr)
		} else if route[routeIdx].Op == trace[i][j] {
			route[routeIdx].RunLen += 1
		} else {
			curr = ByteCigar{RunLen: 1, Op: trace[i][j]}
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
	ReverseBytesCigar(route)
	return m[maxI][maxJ], route, 0, maxI, 0, maxJ
}

func LeftDynamicAln(alpha []dna.Base, beta []dna.Base, scores [][]int64, gapPen int64, m [][]int64, trace [][]CigarOp) (int64, []ByteCigar, int, int, int, int) {
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
			m[i][j], trace[i][j] = ByteMatrixTrace(m[i-1][j-1]+scores[alpha[i-1]][beta[j-1]], m[i][j-1]+gapPen, m[i-1][j]+gapPen)
			if m[i][j] < 0 {
				m[i][j] = 0
			}
		}
	}
	var minI, minJ = len(alpha), len(beta)
	var route []ByteCigar = make([]ByteCigar, 0, 1)
	curr := ByteCigar{}
	//traceback starts in top corner
	for i, j, routeIdx = len(alpha), len(beta), 0; m[i][j] > 0; {
		//if route[routeIdx].RunLength == 0 {
		if len(route) == 0 {
			curr = ByteCigar{RunLen: 1, Op: trace[i][j]}
			route = append(route, curr)
		} else if route[routeIdx].Op == trace[i][j] {
			route[routeIdx].RunLen += 1
		} else {
			curr = ByteCigar{RunLen: 1, Op: trace[i][j]}
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
	ReverseBytesCigar(route)
	return m[len(alpha)][len(beta)], route, minI, len(alpha), minJ, len(beta)
}


func RightAlignTraversal(n *Node, seq []dna.Base, start int, currentPath []uint32, ext int, read []dna.Base, m [][]int64, trace [][]simpleio.CigarOp) ([]simpleio.ByteCigar, int64, int, int, []uint32) {
	path := make([]uint32, 0, len(currentPath))
	copy(path, currentPath)
	var bestTargetEnd, bestQueryEnd, targetEnd, queryEnd int
	var bestScore, score int64
	var bestAlignment, alignment []simpleio.ByteCigar
	var bestPath []uint32

	if len(seq) >= ext {
		log.Fatalf("Error: the length of DNA sequence in previous nodes should not be enough to satisfy the desired extenion.\n")
	}
	var availableBases int = len(seq) + len(n.Seq) - start
	var targetLength int = common.Min(availableBases, ext)
	var basesToTake int = targetLength - len(seq)
	var s []dna.Base = make([]dna.Base, targetLength)
	copy(s[0:len(seq)], seq)
	copy(s[len(seq):targetLength], n.Seq[start:start+basesToTake])

	if availableBases >= ext || len(n.Next) == 0 {
		score, alignment, _, targetEnd, _, queryEnd = simpleio.RightDynamicAln(s, read, HumanChimpTwoScoreMatrix, -600, m, trace)
		return alignment, score, targetEnd+start, queryEnd, path
	} else {
		bestScore = -1
		for _, i := range n.Next {
			alignment, score, targetEnd, queryEnd, path = RightAlignTraversal(i.Dest, s, 0, path, ext, read, m, trace)
			if score > bestScore {
				bestScore = score
				bestAlignment = alignment
				bestTargetEnd = targetEnd
				bestQueryEnd = queryEnd
				bestPath = path
			}
		}
	}
	return bestAlignment, bestScore, bestTargetEnd+start, bestQueryEnd, bestPath
}
// Get Target sequence to align to
func getLeftTargetBases(n *Node, extenion int, tEndPos int, headSeq []dna.Base) []dna.Base {
	var targetLen int = common.Min(len(headSeq)+tEndPos, extenion)
	var basesToTake int = targetLen - len(headSeq)
	ans := make([]dna.Base, targetLen)
	copy(ans[0:basesToTake], n.Seq[tEndPos-basesToTake:tEndPos])
	copy(ans[basesToTake:targetLen], headSeq)
	return ans
}
func LeftAlignTraversal(n *Node, seq []dna.Base, refEnd int, currentPath []uint32, ext int, read []dna.Base, m [][]int64, trace [][]simpleio.CigarOp) ([]simpleio.ByteCigar, int64, int, int, []uint32) {
	path := make([]uint32, 0, len(currentPath))
	copy(path, currentPath)
	AddPath(path, n.Id)
	var bestQueryStart, queryStart, refStart, bestRefStart int
	var bestScore, score int64
	var bestAlignment, alignment []simpleio.ByteCigar
	var bestPath []uint32
	s := getLeftTargetBases(n, ext, refEnd, seq)

	if len(s)-len(seq) >= ext || len(n.Next) == 0 {
		score, alignment, refStart, _, queryStart, _ = simpleio.LeftDynamicAln(s, read, HumanChimpTwoScoreMatrix, -600, m, trace)
		simpleio.ReversePath(path)
		return alignment, score, refEnd - len(s)-len(seq) + refStart, queryStart, path
	} else {
		bestScore = -1

		for _, i := range n.Prev {
			alignment, score, refStart, queryStart, path = LeftAlignTraversal(i.Dest, s, len(i.Dest.Seq), path, ext, read, m, trace)
			if score > bestScore {
				bestScore = score
				bestAlignment = alignment
				bestRefStart = refStart
				bestQueryStart = queryStart
				bestPath = path
			}
		}
	}
	return bestAlignment, bestScore, refEnd - len(s)-len(seq) +bestRefStart, bestQueryStart, bestPath
}

func SoftClipBases(front int, lengthOfRead int, cig []ByteCigar) []ByteCigar {
	var runLen int = QueryRunLen(cig)
	if runLen < lengthOfRead {
		answer := make([]ByteCigar, 0, len(cig)+2)
		if front > 0 {
			answer = append(answer, ByteCigar{RunLen: uint32(front), Op: 'S'})
		}
		answer = append(answer, cig...)
		if front+QueryRunLen(cig) < lengthOfRead {
			answer = append(answer, ByteCigar{RunLen: uint32(lengthOfRead-front - runLen), Op: 'S'})
		}
		return answer
	} else {
		return cig
	}
}

