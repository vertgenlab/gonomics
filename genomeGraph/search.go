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
	route    []cigar.Cigar
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
	leftAlignment  []cigar.Cigar
	rightAlignment []cigar.Cigar
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

func LeftAlignTraversal(n *Node, seq []dna.Base, start int, currentPath []uint32, query []dna.Base, config *GraphSettings, matrix *MatrixAln, sk scoreKeeper, res dynamicScoreKeeper, pool *sync.Pool) ([]cigar.Cigar, int64, int, int, []uint32) {
	cache := pool.Get().(*dnaPool)
	defer pool.Put(cache)
	cache.Seq, cache.Path = cache.Seq[:0], cache.Path[:0]
	cache.Seq = getLeftTargetBases(n, config.Extention, start, seq, cache.Seq)
	cache.Path = make([]uint32, len(currentPath))

	copy(cache.Path, currentPath)
	AddPath(cache.Path, n.Id)

	currLen := len(seq)
	if currLen+start >= config.Extention || len(n.Prev) == 0 {
		sk.leftScore, sk.leftAlignment, sk.targetStart, sk.queryStart = LeftDynamicAln(cache.Seq, query, config.ScoreMatrix, matrix, -600, res)
		sk.targetStart = start - len(cache.Seq) - currLen + sk.targetStart
		sk.leftPath = cache.Path
		return sk.leftAlignment, sk.leftScore, sk.targetStart, sk.queryStart, sk.leftPath
	} else {
		//A very negative number
		sk.leftScore = math.MinInt64
		for _, i := range n.Prev {
			res.route, cache.currScore, cache.targetStart, cache.queryStart, cache.Path = LeftAlignTraversal(i.Dest, cache.Seq, len(i.Dest.Seq), cache.Path, query, config, matrix, sk, res, pool)
			if cache.currScore > sk.leftScore {
				sk.leftScore = cache.currScore
				sk.leftAlignment = res.route
				sk.targetStart = start - len(cache.Seq) - currLen + cache.targetStart
				sk.queryStart = cache.queryStart
				sk.leftPath = cache.Path
			}
		}

		cigar.ReverseCigar(sk.leftAlignment)
		reversePath(sk.leftPath)
		return sk.leftAlignment, sk.leftScore, sk.targetStart, sk.queryStart, sk.leftPath
	}
}

func RightAlignTraversal(n *Node, seq []dna.Base, end int, currentPath []uint32, query []dna.Base, config *GraphSettings, matrix *MatrixAln, sk scoreKeeper, res dynamicScoreKeeper, pool *sync.Pool) ([]cigar.Cigar, int64, int, int, []uint32) {
	cache := pool.Get().(*dnaPool)
	defer pool.Put(cache)

	cache.Seq, cache.Path = cache.Seq[:0], cache.Path[:0]
	cache.Seq = getRightBases(n, config.Extention, end, seq, cache.Seq)
	cache.Path = make([]uint32, len(currentPath))
	copy(cache.Path, currentPath)
	if len(seq)+len(n.Seq)-end >= config.Extention || len(n.Next) == 0 {
		sk.rightScore, sk.rightAlignment, sk.targetEnd, sk.queryEnd = RightDynamicAln(cache.Seq, query, config, matrix, res)
		sk.rightPath = cache.Path
		return sk.rightAlignment, sk.rightScore, sk.targetEnd + end, sk.queryEnd, sk.rightPath
	} else {
		sk.rightScore = math.MinInt64
		for _, i := range n.Next {
			res.route, cache.currScore, cache.targetEnd, cache.queryEnd, cache.Path = RightAlignTraversal(i.Dest, cache.Seq, 0, cache.Path, query, config, matrix, sk, res, pool)
			if cache.currScore > sk.rightScore {
				sk.rightScore = cache.currScore
				sk.rightAlignment = res.route
				sk.targetEnd = cache.targetEnd
				sk.queryEnd = cache.queryEnd
				sk.rightPath = cache.Path
			}
		}
		cigar.ReverseCigar(sk.rightAlignment)
		return sk.rightAlignment, sk.rightScore, sk.targetEnd + end, sk.queryEnd, sk.rightPath
	}
}

func LeftDynamicAln(alpha []dna.Base, beta []dna.Base, scores [][]int64, matrix *MatrixAln, gapPen int64, res dynamicScoreKeeper) (int64, []cigar.Cigar, int, int) {
	resetDynamicScore(res)
	columns, rows := len(alpha), len(beta)

	for res.i = 0; res.i <= columns; res.i++ {
		matrix.m[res.i][0] = 0
	}

	for res.j = 0; res.j <= rows; res.j++ {
		matrix.m[0][res.j] = 0
	}

	for res.i = 1; res.i <= columns; res.i++ {
		for res.j = 1; res.j <= rows; res.j++ {
			matrix.m[res.i][res.j], matrix.trace[res.i][res.j] = cigar.TripleMaxTrace(matrix.m[res.i-1][res.j-1]+scores[alpha[res.i-1]][beta[res.j-1]], matrix.m[res.i][res.j-1]+gapPen, matrix.m[res.i-1][res.j]+gapPen)
			matrix.m[res.i][res.j] = numbers.Max(matrix.m[res.i][res.j], 0)
		}
	}

	for res.i, res.j, res.routeIdx = columns, rows, 0; matrix.m[res.i][res.j] > 0; {
		if len(res.route) == 0 || res.route[len(res.route)-1].Op != matrix.trace[res.i][res.j] {
			res.route = append(res.route, cigar.Cigar{RunLength: 1, Op: matrix.trace[res.i][res.j]})
		} else {
			res.route[len(res.route)-1].RunLength++
		}
		switch matrix.trace[res.i][res.j] {
		case cigar.Match:
			res.i, res.j = res.i-1, res.j-1
		case cigar.Insertion:
			res.j -= 1
		case cigar.Deletion:
			res.i -= 1
		default:
			log.Fatalf("Error: unexpected traceback %c\n", matrix.trace[res.i][res.j])
		}
	}
	return matrix.m[columns][rows], res.route, res.i, res.j
}

func RightDynamicAln(alpha []dna.Base, beta []dna.Base, config *GraphSettings, matrix *MatrixAln, res dynamicScoreKeeper) (int64, []cigar.Cigar, int, int) {
	columns, rows := len(alpha), len(beta)
	resetDynamicScore(res)
	var maxI int
	var maxJ int

	for res.i = 0; res.i <= columns; res.i++ {
		matrix.m[res.i][0] = config.GapPenalty * int64(res.i)
		matrix.trace[res.i][0] = cigar.Deletion // Initialize first column traces

	}
	for res.j = 0; res.j <= rows; res.j++ {
		matrix.m[0][res.j] = config.GapPenalty * int64(res.j)
		matrix.trace[0][res.j] = cigar.Insertion // Initialize first row traces
	}

	for res.i = 1; res.i <= columns; res.i++ {
		for res.j = 1; res.j <= rows; res.j++ {
			matrix.m[res.i][res.j], matrix.trace[res.i][res.j] = cigar.TripleMaxTrace(
				matrix.m[res.i-1][res.j-1]+config.ScoreMatrix[alpha[res.i-1]][beta[res.j-1]],
				matrix.m[res.i][res.j-1]+config.GapPenalty,
				matrix.m[res.i-1][res.j]+config.GapPenalty,
			)
			if matrix.m[res.i][res.j] > res.currMax {
				res.currMax = matrix.m[res.i][res.j]
				maxI = res.i
				maxJ = res.j
			}
		}
	}
	for res.i, res.j, res.routeIdx = maxI, maxJ, 0; res.i > 0 || res.j > 0; {
		if len(res.route) == 0 || res.route[len(res.route)-1].Op != matrix.trace[res.i][res.j] {
			res.route = append(res.route, cigar.Cigar{RunLength: 1, Op: matrix.trace[res.i][res.j]})
		} else {
			res.route[len(res.route)-1].RunLength++
		}
		switch matrix.trace[res.i][res.j] {
		case cigar.Match:
			res.i, res.j = res.i-1, res.j-1
		case cigar.Insertion:
			res.j -= 1
		case cigar.Deletion:
			res.i -= 1
		default:
			log.Fatalf("Error: unexpected traceback with %c\n", matrix.trace[res.i][res.j])
		}
	}
	return matrix.m[maxI][maxJ], res.route, maxI, maxJ
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
