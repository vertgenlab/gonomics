package genomeGraph

import (
	"bytes"
	"io"
	"log"
	"math"
	"sort"
	"sync"

	"github.com/vertgenlab/gonomics/cigar"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/giraf"
	"github.com/vertgenlab/gonomics/numbers"
)

const (
	defaultMatrixSize int = 2480
)

type memoryPool struct {
	Hits   []Seed
	Worker []Seed
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
	tailSeed     Seed

	currSeed       Seed
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
				Hits:   make([]Seed, 0, 10000),
				Worker: make([]Seed, 0, 10000),
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
	currLen := len(seq)
	return append(append(ans, n.Seq[refEnd-numbers.Min(currLen+refEnd, extension)-currLen:refEnd]...), seq...)
}

func getRightBases(n *Node, extension int, start int, seq []dna.Base, ans []dna.Base) []dna.Base {
	currLen := len(seq)
	return append(append(ans, seq...), n.Seq[start:start+numbers.Min(currLen+len(n.Seq)-start, extension)-len(seq)]...)
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

func LeftAlignTraversal(n *Node, seq []dna.Base, start int, currentPath []uint32, query []dna.Base, config *GraphSettings, matrix *MatrixAln, sk scoreKeeper, dynamicScore dynamicScoreKeeper, pool *sync.Pool) ([]cigar.ByteCigar, int64, int, int, []uint32) {
	s := pool.Get().(*dnaPool)
	s.Seq, s.Path = s.Seq[:0], s.Path[:0]
	s.Seq = getLeftTargetBases(n, config.Extention, start, seq, s.Seq)
	s.Path = make([]uint32, len(currentPath))

	copy(s.Path, currentPath)
	AddPath(s.Path, n.Id)

	currLen := len(seq)
	if currLen+start >= config.Extention || len(n.Prev) == 0 {
		sk.leftScore, sk.leftAlignment, sk.targetStart, sk.queryStart = LeftDynamicAln(s.Seq, query, config.ScoreMatrix, matrix, -600, dynamicScore)
		sk.targetStart = start - len(s.Seq) - currLen + sk.targetStart
		sk.leftPath = s.Path
		pool.Put(s)
		return sk.leftAlignment, sk.leftScore, sk.targetStart, sk.queryStart, sk.leftPath
	} else {
		//A very negative number
		sk.leftScore = math.MinInt64
		for _, i := range n.Prev {
			dynamicScore.route, s.currScore, s.targetStart, s.queryStart, s.Path = LeftAlignTraversal(i.Dest, s.Seq, len(i.Dest.Seq), s.Path, query, config, matrix, sk, dynamicScore, pool)
			if s.currScore > sk.leftScore {
				sk.leftScore = s.currScore
				sk.leftAlignment = dynamicScore.route
				sk.targetStart = start - len(s.Seq) - currLen + s.targetStart
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

func RightAlignTraversal(n *Node, seq []dna.Base, end int, currentPath []uint32, query []dna.Base, config *GraphSettings, matrix *MatrixAln, sk scoreKeeper, dynamicScore dynamicScoreKeeper, pool *sync.Pool) ([]cigar.ByteCigar, int64, int, int, []uint32) {
	s := pool.Get().(*dnaPool)
	s.Seq, s.Path = s.Seq[:0], s.Path[:0]
	s.Seq = getRightBases(n, config.Extention, end, seq, s.Seq)
	s.Path = make([]uint32, len(currentPath))
	copy(s.Path, currentPath)
	if len(seq)+len(n.Seq)-end >= config.Extention || len(n.Next) == 0 {
		sk.rightScore, sk.rightAlignment, sk.targetEnd, sk.queryEnd = RightDynamicAln(s.Seq, query, config, matrix, dynamicScore)
		sk.rightPath = s.Path
		pool.Put(s)
		return sk.rightAlignment, sk.rightScore, sk.targetEnd + end, sk.queryEnd, sk.rightPath
	} else {
		sk.rightScore = math.MinInt64
		for _, i := range n.Next {
			dynamicScore.route, s.currScore, s.targetEnd, s.queryEnd, s.Path = RightAlignTraversal(i.Dest, s.Seq, 0, s.Path, query, config, matrix, sk, dynamicScore, pool)
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
		return sk.rightAlignment, sk.rightScore, sk.targetEnd + end, sk.queryEnd, sk.rightPath
	}
}

func LeftDynamicAln(alpha []dna.Base, beta []dna.Base, scores [][]int64, matrix *MatrixAln, gapPen int64, dynamicScore dynamicScoreKeeper) (int64, []cigar.ByteCigar, int, int) {
	resetDynamicScore(dynamicScore)
	columns, rows := len(alpha), len(beta)
	for dynamicScore.i = 0; dynamicScore.i < columns+1; dynamicScore.i++ {
		matrix.m[dynamicScore.i][0] = 0
	}

	for dynamicScore.j = 0; dynamicScore.j < rows+1; dynamicScore.j++ {
		matrix.m[0][dynamicScore.j] = 0
	}

	for dynamicScore.i = 1; dynamicScore.i < columns+1; dynamicScore.i++ {
		for dynamicScore.j = 1; dynamicScore.j < rows+1; dynamicScore.j++ {
			matrix.m[dynamicScore.i][dynamicScore.j], matrix.trace[dynamicScore.i][dynamicScore.j] = cigar.ByteMatrixTrace(matrix.m[dynamicScore.i-1][dynamicScore.j-1]+scores[alpha[dynamicScore.i-1]][beta[dynamicScore.j-1]], matrix.m[dynamicScore.i][dynamicScore.j-1]+gapPen, matrix.m[dynamicScore.i-1][dynamicScore.j]+gapPen)
			if matrix.m[dynamicScore.i][dynamicScore.j] < 0 {
				matrix.m[dynamicScore.i][dynamicScore.j] = 0
			}
		}
	}

	for dynamicScore.i, dynamicScore.j, dynamicScore.routeIdx = columns, rows, 0; matrix.m[dynamicScore.i][dynamicScore.j] > 0; {
		if len(dynamicScore.route) == 0 {
			dynamicScore.route = append(dynamicScore.route, cigar.ByteCigar{RunLen: 1, Op: matrix.trace[dynamicScore.i][dynamicScore.j]})
		} else if dynamicScore.route[dynamicScore.routeIdx].Op == matrix.trace[dynamicScore.i][dynamicScore.j] {
			dynamicScore.route[dynamicScore.routeIdx].RunLen += 1
		} else {
			dynamicScore.route = append(dynamicScore.route, cigar.ByteCigar{RunLen: 1, Op: matrix.trace[dynamicScore.i][dynamicScore.j]})
			dynamicScore.routeIdx++
		}
		switch matrix.trace[dynamicScore.i][dynamicScore.j] {
		case cigar.Match:
			dynamicScore.i, dynamicScore.j = dynamicScore.i-1, dynamicScore.j-1
		case cigar.Insertion:
			dynamicScore.j -= 1
		case cigar.Deletion:
			dynamicScore.i -= 1
		default:
			log.Fatalf("Error: unexpected traceback %c\n", matrix.trace[dynamicScore.i][dynamicScore.j])
		}
	}
	return matrix.m[columns][rows], dynamicScore.route, dynamicScore.i, dynamicScore.j
}

func RightDynamicAln(alpha []dna.Base, beta []dna.Base, config *GraphSettings, matrix *MatrixAln, dynamicScore dynamicScoreKeeper) (int64, []cigar.ByteCigar, int, int) {
	columns, rows := len(alpha)+1, len(beta)+1
	resetDynamicScore(dynamicScore)
	var maxI int
	var maxJ int
	for dynamicScore.i = 0; dynamicScore.i < columns; dynamicScore.i++ {
		for dynamicScore.j = 0; dynamicScore.j < rows; dynamicScore.j++ {
			if dynamicScore.i == 0 && dynamicScore.j == 0 {
				matrix.m[dynamicScore.i][dynamicScore.j] = 0
			} else if dynamicScore.i == 0 {
				matrix.m[dynamicScore.i][dynamicScore.j] = matrix.m[dynamicScore.i][dynamicScore.j-1] + config.GapPenalty
				matrix.trace[dynamicScore.i][dynamicScore.j] = cigar.Insertion
			} else if dynamicScore.j == 0 {
				matrix.m[dynamicScore.i][dynamicScore.j] = matrix.m[dynamicScore.i-1][dynamicScore.j] + config.GapPenalty
				matrix.trace[dynamicScore.i][dynamicScore.j] = cigar.Deletion
			} else {
				matrix.m[dynamicScore.i][dynamicScore.j], matrix.trace[dynamicScore.i][dynamicScore.j] = cigar.ByteMatrixTrace(matrix.m[dynamicScore.i-1][dynamicScore.j-1]+config.ScoreMatrix[alpha[dynamicScore.i-1]][beta[dynamicScore.j-1]], matrix.m[dynamicScore.i][dynamicScore.j-1]+config.GapPenalty, matrix.m[dynamicScore.i-1][dynamicScore.j]+config.GapPenalty)
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
		case cigar.Match:
			dynamicScore.i, dynamicScore.j = dynamicScore.i-1, dynamicScore.j-1
		case cigar.Insertion:
			dynamicScore.j -= 1
		case cigar.Deletion:
			dynamicScore.i -= 1
		default:
			log.Fatalf("Error: unexpected traceback with %c\n", matrix.trace[dynamicScore.i][dynamicScore.j])
		}
	}
	return matrix.m[maxI][maxJ], dynamicScore.route, maxI, maxJ
}

func ReversePath(alpha []uint32) {
	var i, off int
	for i = len(alpha)/2 - 1; i >= 0; i-- {
		off = len(alpha) - 1 - i
		alpha[i], alpha[off] = alpha[off], alpha[i]
	}
}

type seedHelper struct {
	currHits                                  []uint64
	codedNodeCoord                            uint64
	seqKey                                    uint64
	keyShift                                  uint
	keyIdx, keyOffset, readOffset, nodeOffset int
	nodeIdx, nodePos                          int64
	leftMatches                               int
	rightMatches                              int
	tempSeed                                  Seed
}

func newSeedBuilder() *seedHelper {
	var tmp Seed = Seed{}
	return &seedHelper{
		currHits: make([]uint64, 0, 20),
		tempSeed: tmp,
	}
}

func restartSeedHelper(helper *seedHelper) {
	helper.currHits = helper.currHits[:0]
	helper.keyIdx, helper.keyOffset, helper.readOffset, helper.nodeOffset = 0, 0, 0, 0
	helper.nodeIdx, helper.nodePos = 0, 0
	helper.seqKey, helper.codedNodeCoord = 0, 0
	helper.leftMatches = 0
}

func SortSeedLen(seeds []Seed) {
	sort.Slice(seeds, func(i, j int) bool { return seeds[i].TotalLength > seeds[j].TotalLength })
}

func SoftClipBases(front int, lengthOfRead int, cig []cigar.ByteCigar) []cigar.ByteCigar {
	var runLen int = cigar.QueryRunLen(cig)
	if runLen < lengthOfRead {
		answer := make([]cigar.ByteCigar, 0, len(cig)+2)
		if front > 0 {
			answer = append(answer, cigar.ByteCigar{RunLen: uint16(front), Op: 'S'})
		}
		answer = append(answer, cig...)
		if front+cigar.QueryRunLen(cig) < lengthOfRead {
			answer = append(answer, cigar.ByteCigar{RunLen: uint16(lengthOfRead - front - runLen), Op: 'S'})
		}
		return answer
	} else {
		return cig
	}
}

func SimpleWriteGirafPair(filename string, input <-chan giraf.GirafPair, wg *sync.WaitGroup) {
	file := fileio.EasyCreate(filename)
	var buf *bytes.Buffer
	var simplePool = sync.Pool{
		New: func() interface{} {
			return &bytes.Buffer{}
		},
	}
	var err error
	for gp := range input {
		buf = simplePool.Get().(*bytes.Buffer)
		buf.Reset()
		_, err = buf.WriteString(giraf.ToString(&gp.Fwd))
		exception.PanicOnErr(err)
		err = buf.WriteByte('\n')
		exception.PanicOnErr(err)
		_, err = buf.WriteString(giraf.ToString(&gp.Rev))
		exception.PanicOnErr(err)
		err = buf.WriteByte('\n')
		exception.PanicOnErr(err)
		_, err = io.Copy(file, buf)
		exception.PanicOnErr(err)
		simplePool.Put(buf)
	}
	err = file.Close()
	exception.PanicOnErr(err)
	wg.Done()
}
