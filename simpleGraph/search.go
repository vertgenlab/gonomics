package simpleGraph

import (
	"github.com/vertgenlab/gonomics/cigar"
	"github.com/vertgenlab/gonomics/common"
	"github.com/vertgenlab/gonomics/dna"
	"log"
	"math"
	"sync"
)

const (
	defaultMatrixSize int = 10000
)

type memoryPool struct {
	Hits   []*SeedDev
	Worker []*SeedDev
}

type MatrixAln struct {
	m     [][]int64
	trace [][]byte
}

type dnaPool struct {
	Seq         []dna.Base
	Path        []uint32
	queryStart  int
	queryEnd    int
	targetStart int
	targetEnd   int
	currScore   int64
}

type dynamicScoreKeeper struct {
	i        int
	j        int
	routeIdx int
	currMax  int64
	route    []cigar.ByteCigar
}

type scoreKeeper struct {
	targetStart  int
	targetEnd    int
	queryStart   int
	queryEnd     int
	extension    int
	currScore    int64
	seedScore    int64
	perfectScore int64
	leftScore    int64
	rightScore   int64
	leftPath     []uint32
	rightPath    []uint32
	leftSeq      []dna.Base
	rightSeq     []dna.Base
	currSeq      []dna.Base
	tailSeed     *SeedDev

	currSeed       *SeedDev
	leftAlignment  []cigar.ByteCigar
	rightAlignment []cigar.ByteCigar
}

func NewDnaPool() sync.Pool {
	return sync.Pool{
		New: func() interface{} {
			dnaSeq := dnaPool{
				Seq:         make([]dna.Base, 0, 150),
				Path:        make([]uint32, 0, 10),
				queryStart:  0,
				targetStart: 0,
				targetEnd:   0,
				queryEnd:    0,
			}
			return &dnaSeq
		},
	}
}

func NewMemSeedPool() sync.Pool {
	return sync.Pool{
		New: func() interface{} {
			pool := memoryPool{
				Hits:   make([]*SeedDev, 10000),
				Worker: make([]*SeedDev, 10000),
			}
			return &pool
		},
	}
}

func resetDynamicScore(sk dynamicScoreKeeper) {
	sk.route = sk.route[:0]
	sk.currMax = 0
}
func NewSwMatrix(size int) MatrixAln {
	sw := MatrixAln{}
	sw.m, sw.trace = MatrixSetup(size)
	return sw
}

func MatrixSetup(size int) ([][]int64, [][]byte) {
	m := make([][]int64, size)
	trace := make([][]byte, size)
	for idx := range m {
		m[idx] = make([]int64, size)
		trace[idx] = make([]byte, size)
	}
	return m, trace
}

func resetScoreKeeper(sk scoreKeeper) {
	sk.targetStart, sk.targetEnd = 0, 0
	sk.queryStart, sk.queryEnd = 0, 0
	sk.currScore, sk.seedScore = 0, 0
	sk.perfectScore = 0
	sk.leftAlignment, sk.rightAlignment = sk.leftAlignment[:0], sk.rightAlignment[:0]
	sk.leftPath, sk.rightPath = sk.leftPath[:0], sk.rightPath[:0]
	sk.leftSeq, sk.rightSeq, sk.currSeq = sk.leftSeq[:0], sk.rightSeq[:0], sk.currSeq[:0]
	sk.leftScore, sk.rightScore = 0, 0
}

func getLeftTargetBases(n *Node, extension int, refEnd int, seq []dna.Base, ans []dna.Base) []dna.Base {
	var availableBases int = len(seq) + refEnd
	var targetLength int = common.Min(availableBases, extension)
	var basesToTake int = targetLength - len(seq)
	return append(append(ans, n.Seq[refEnd-basesToTake:refEnd]...), seq...)
}

func getRightBases(n *Node, extension int, start int, seq []dna.Base, ans []dna.Base) []dna.Base {
	var availableBases int = len(seq) + len(n.Seq) - start
	var targetLength int = common.Min(availableBases, extension)
	var basesToTake int = targetLength - len(seq)
	return append(append(ans, seq...), n.Seq[start:start+basesToTake]...)
}

/*
func leftBasesFromTwoBit(n *Node, extension int, refEnd int, seq []dna.Base, ans []dna.Base) []dna.Base {
	var availableBases int = len(seq) + refEnd
	var targetLength int = common.Min(availableBases, extension)
	var basesToTake int = targetLength - len(seq)
	return append(append(ans, seq...), dnaTwoBit.GetFrag(n.SeqTwoBit, refEnd-basesToTake, refEnd)...)
}*/

/*
func rightBasesFromTwoBit(n *Node, extension int, start int, seq []dna.Base, ans []dna.Base) []dna.Base {
	var availableBases int = len(seq) + len(n.Seq) - start
	var targetLength int = common.Min(availableBases, extension)
	var basesToTake int = targetLength - len(seq)

	return append(append(ans, seq...), dnaTwoBit.GetFrag(n.SeqTwoBit, start, start+basesToTake)...)
}*/

func LeftAlignTraversal(n *Node, seq []dna.Base, refEnd int, currentPath []uint32, extension int, read []dna.Base, scores [][]int64, matrix *MatrixAln, sk scoreKeeper, dynamicScore dynamicScoreKeeper, pool *sync.Pool) ([]cigar.ByteCigar, int64, int, int, []uint32) {
	if len(seq) >= extension {
		log.Fatalf("Error: left traversal, the length=%d of DNA sequence in previous nodes should not be enough to satisfy the desired extenion=%d.\n", len(seq), extension)
	}
	s := pool.Get().(*dnaPool)
	s.Seq, s.Path = s.Seq[:0], s.Path[:0]
	s.Seq = getLeftTargetBases(n, extension, refEnd, seq, s.Seq)
	s.Path = make([]uint32, len(currentPath))
	copy(s.Path, currentPath)
	AddPath(s.Path, n.Id)
	if len(seq)+refEnd >= extension || len(n.Prev) == 0 {
		sk.leftScore, sk.leftAlignment, sk.targetStart, sk.queryStart = LeftDynamicAln(s.Seq, read, scores, matrix, -600, dynamicScore)
		sk.targetStart = refEnd - len(s.Seq) - len(seq) + sk.targetStart
		sk.leftPath = s.Path
		pool.Put(s)
		return sk.leftAlignment, sk.leftScore, sk.targetStart, sk.queryStart, sk.leftPath
	} else {
		//A very negative number
		sk.leftScore = math.MinInt64
		for _, i := range n.Prev {
			dynamicScore.route, s.currScore, s.targetStart, s.queryStart, s.Path = LeftAlignTraversal(i.Dest, s.Seq, len(i.Dest.Seq), s.Path, extension, read, scores, matrix, sk, dynamicScore, pool)
			if s.currScore > sk.leftScore {
				sk.leftScore = s.currScore
				sk.leftAlignment = dynamicScore.route
				sk.targetStart = refEnd - len(s.Seq) - len(seq) + s.targetStart
				sk.queryStart = s.queryStart
				sk.leftPath = s.Path
			}
		}
		pool.Put(s)
		cigar.ReverseBytesCigar(sk.leftAlignment)
		ReversePath(sk.leftPath)
		return sk.leftAlignment, sk.leftScore, sk.targetStart, sk.queryStart, sk.leftPath
	}
}

func RightAlignTraversal(n *Node, seq []dna.Base, start int, currentPath []uint32, extension int, read []dna.Base, scoreMatrix [][]int64, matrix *MatrixAln, sk scoreKeeper, dynamicScore dynamicScoreKeeper, pool *sync.Pool) ([]cigar.ByteCigar, int64, int, int, []uint32) {
	if len(seq) >= extension {
		log.Fatalf("Error: right traversal, the length=%d of DNA sequence in previous nodes should not be enough to satisfy the desired extenion=%d.\n", len(seq), extension)
	}
	s := pool.Get().(*dnaPool)
	s.Seq, s.Path = s.Seq[:0], s.Path[:0]
	s.Seq = getRightBases(n, extension, start, seq, s.Seq)
	s.Path = make([]uint32, len(currentPath))
	copy(s.Path, currentPath)
	if len(seq)+len(n.Seq)-start >= extension || len(n.Next) == 0 {
		sk.rightScore, sk.rightAlignment, sk.targetEnd, sk.queryEnd = RightDynamicAln(s.Seq, read, scoreMatrix, matrix, -600, dynamicScore)
		sk.rightPath = s.Path
		pool.Put(s)
		return sk.rightAlignment, sk.rightScore, sk.targetEnd + start, sk.queryEnd, sk.rightPath
	} else {
		sk.rightScore = math.MinInt64
		for _, i := range n.Next {
			dynamicScore.route, s.currScore, s.targetEnd, s.queryEnd, s.Path = RightAlignTraversal(i.Dest, s.Seq, 0, s.Path, extension, read, scoreMatrix, matrix, sk, dynamicScore, pool)
			if s.currScore > sk.rightScore {
				sk.rightScore = s.currScore
				sk.rightAlignment = dynamicScore.route
				sk.targetEnd = s.targetEnd
				sk.queryEnd = s.queryEnd
				sk.rightPath = s.Path
			}
		}
		pool.Put(s)
		cigar.ReverseBytesCigar(sk.rightAlignment)
		return sk.rightAlignment, sk.rightScore, sk.targetEnd + start, sk.queryEnd, sk.rightPath
	}
}

func LeftDynamicAln(alpha []dna.Base, beta []dna.Base, scores [][]int64, matrix *MatrixAln, gapPen int64, dynamicScore dynamicScoreKeeper) (int64, []cigar.ByteCigar, int, int) {
	resetDynamicScore(dynamicScore)
	for dynamicScore.i = 0; dynamicScore.i < len(alpha)+1; dynamicScore.i++ {
		matrix.m[dynamicScore.i][0] = 0
	}

	for dynamicScore.j = 0; dynamicScore.j < len(beta)+1; dynamicScore.j++ {
		matrix.m[0][dynamicScore.j] = 0
	}

	for dynamicScore.i = 1; dynamicScore.i < len(alpha)+1; dynamicScore.i++ {
		for dynamicScore.j = 1; dynamicScore.j < len(beta)+1; dynamicScore.j++ {
			matrix.m[dynamicScore.i][dynamicScore.j], matrix.trace[dynamicScore.i][dynamicScore.j] = cigar.ByteMatrixTrace(matrix.m[dynamicScore.i-1][dynamicScore.j-1]+scores[alpha[dynamicScore.i-1]][beta[dynamicScore.j-1]], matrix.m[dynamicScore.i][dynamicScore.j-1]+gapPen, matrix.m[dynamicScore.i-1][dynamicScore.j]+gapPen)
			if matrix.m[dynamicScore.i][dynamicScore.j] < 0 {
				matrix.m[dynamicScore.i][dynamicScore.j] = 0
			}
		}
	}

	for dynamicScore.i, dynamicScore.j, dynamicScore.routeIdx = len(alpha), len(beta), 0; matrix.m[dynamicScore.i][dynamicScore.j] > 0; {
		if len(dynamicScore.route) == 0 {
			dynamicScore.route = append(dynamicScore.route, cigar.ByteCigar{RunLen: 1, Op: matrix.trace[dynamicScore.i][dynamicScore.j]})
		} else if dynamicScore.route[dynamicScore.routeIdx].Op == matrix.trace[dynamicScore.i][dynamicScore.j] {
			dynamicScore.route[dynamicScore.routeIdx].RunLen += 1
		} else {
			dynamicScore.route = append(dynamicScore.route, cigar.ByteCigar{RunLen: 1, Op: matrix.trace[dynamicScore.i][dynamicScore.j]})
			dynamicScore.routeIdx++
		}
		switch matrix.trace[dynamicScore.i][dynamicScore.j] {
		case 'M':
			dynamicScore.i, dynamicScore.j = dynamicScore.i-1, dynamicScore.j-1
		case 'I':
			dynamicScore.j -= 1
		case 'D':
			dynamicScore.i -= 1
		default:
			log.Fatalf("Error: unexpected traceback %c\n", matrix.trace[dynamicScore.i][dynamicScore.j])
		}

	}
	return matrix.m[len(alpha)][len(beta)], dynamicScore.route, dynamicScore.i, dynamicScore.j
}

func RightDynamicAln(alpha []dna.Base, beta []dna.Base, scores [][]int64, matrix *MatrixAln, gapPen int64, dynamicScore dynamicScoreKeeper) (int64, []cigar.ByteCigar, int, int) {
	resetDynamicScore(dynamicScore)
	var maxI int
	var maxJ int
	for dynamicScore.i = 0; dynamicScore.i < len(alpha)+1; dynamicScore.i++ {
		for dynamicScore.j = 0; dynamicScore.j < len(beta)+1; dynamicScore.j++ {
			if dynamicScore.i == 0 && dynamicScore.j == 0 {
				matrix.m[dynamicScore.i][dynamicScore.j] = 0
			} else if dynamicScore.i == 0 {
				matrix.m[dynamicScore.i][dynamicScore.j] = matrix.m[dynamicScore.i][dynamicScore.j-1] + gapPen
				matrix.trace[dynamicScore.i][dynamicScore.j] = 'I'
			} else if dynamicScore.j == 0 {
				matrix.m[dynamicScore.i][dynamicScore.j] = matrix.m[dynamicScore.i-1][dynamicScore.j] + gapPen
				matrix.trace[dynamicScore.i][dynamicScore.j] = 'D'
			} else {
				matrix.m[dynamicScore.i][dynamicScore.j], matrix.trace[dynamicScore.i][dynamicScore.j] = cigar.ByteMatrixTrace(matrix.m[dynamicScore.i-1][dynamicScore.j-1]+scores[alpha[dynamicScore.i-1]][beta[dynamicScore.j-1]], matrix.m[dynamicScore.i][dynamicScore.j-1]+gapPen, matrix.m[dynamicScore.i-1][dynamicScore.j]+gapPen)
			}
			if matrix.m[dynamicScore.i][dynamicScore.j] > dynamicScore.currMax {
				dynamicScore.currMax = matrix.m[dynamicScore.i][dynamicScore.j]
				maxI = dynamicScore.i
				maxJ = dynamicScore.j
			}
		}
	}
	for dynamicScore.i, dynamicScore.j, dynamicScore.routeIdx = maxI, maxJ, 0; dynamicScore.i > 0 || dynamicScore.j > 0; {
		if len(dynamicScore.route) == 0 {
			dynamicScore.route = append(dynamicScore.route, cigar.ByteCigar{RunLen: 1, Op: matrix.trace[dynamicScore.i][dynamicScore.j]})
		} else if dynamicScore.route[dynamicScore.routeIdx].Op == matrix.trace[dynamicScore.i][dynamicScore.j] {
			dynamicScore.route[dynamicScore.routeIdx].RunLen += 1
		} else {
			dynamicScore.route = append(dynamicScore.route, cigar.ByteCigar{RunLen: 1, Op: matrix.trace[dynamicScore.i][dynamicScore.j]})
			dynamicScore.routeIdx++
		}
		switch matrix.trace[dynamicScore.i][dynamicScore.j] {
		case 'M':
			dynamicScore.i, dynamicScore.j = dynamicScore.i-1, dynamicScore.j-1
		case 'I':
			dynamicScore.j -= 1
		case 'D':
			dynamicScore.i -= 1
		default:
			log.Fatalf("Error: unexpected traceback with %c\n", matrix.trace[dynamicScore.i][dynamicScore.j])
		}
	}
	return matrix.m[maxI][maxJ], dynamicScore.route, maxI, maxJ
}

func ReversePath(alpha []uint32) {
	var i, off int
	for i, off = len(alpha)/2-1, len(alpha)-1-i; i >= 0; i-- {
		alpha[i], alpha[off] = alpha[off], alpha[i]
	}
}
