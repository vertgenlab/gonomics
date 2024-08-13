package genomeGraph

import (
	"bytes"
	"io"
	"math"
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

func getLeftBases(n *Node, extension int, refEnd int, ans []dna.Base) []dna.Base {
	seqLen := len(ans)
	var basesToTake int = (numbers.Min(seqLen+refEnd, extension) - seqLen)

	ans = append(ans, make([]dna.Base, basesToTake)...)
	copy(ans[seqLen:], n.Seq[refEnd-basesToTake:refEnd])
	return ans
}

func getRightBases(n *Node, extension, refStart int, ans []dna.Base) []dna.Base {
	seqLen := len(ans)
	basesToTake := numbers.Min(len(n.Seq)-refStart, extension-seqLen)

	// Preallocate for efficiency
	ans = append(ans, make([]dna.Base, basesToTake)...)

	// Copy the bases from n.Seq
	copy(ans[seqLen:], n.Seq[refStart:refStart+basesToTake])

	return ans
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


func LeftAlignTraversal(n *Node, seq []dna.Base, refEnd int, currentPath []uint32, read []dna.Base, extension int, config *GraphSettings, matrix *sync.Pool, sk scoreKeeper, dynamicScore dynamicScoreKeeper, pool *sync.Pool) ([]cigar.Cigar, int64, int, int, []uint32) {
	s := pool.Get().(*dnaPool)
	defer pool.Put(s)
	s.Seq, s.Path = s.Seq[:0], s.Path[:0]
	s.Seq = getLeftBases(n, extension, refEnd, seq)
	s.Path = make([]uint32, len(currentPath))
	copy(s.Path, currentPath)
	AddPath(s.Path, n.Id)
	if len(s.Seq)+refEnd >= extension || len(n.Prev) == 0 {
		sk.leftScore, sk.leftAlignment, sk.targetStart, _, sk.queryStart, _ = LeftLocal(s.Seq, read, config, matrix)
		sk.targetStart = refEnd - len(s.Seq) + sk.targetStart
		sk.leftPath = s.Path
		return sk.leftAlignment, sk.leftScore, sk.targetStart, sk.queryStart, sk.leftPath
	} else {
		//A very negative number
		sk.leftScore = math.MinInt64
		for _, i := range n.Prev {
			dynamicScore.route, s.currScore, s.targetStart, s.queryStart, s.Path = LeftAlignTraversal(i.Dest, s.Seq, len(i.Dest.Seq), s.Path, read, extension, config, matrix, sk, dynamicScore, pool)
			if s.currScore > sk.leftScore {
				sk.leftScore = s.currScore
				sk.leftAlignment = dynamicScore.route
				sk.targetStart = refEnd - len(s.Seq) + s.targetStart
				sk.queryStart = s.queryStart
				sk.leftPath = s.Path
			}
		}
		
		cigar.ReverseCigar(sk.leftAlignment)
		ReversePath(sk.leftPath)
		return sk.leftAlignment, sk.leftScore, sk.targetStart, sk.queryStart, sk.leftPath
	}
}

func RightAlignTraversal(n *Node, seq []dna.Base, start int, currentPath []uint32, read []dna.Base, extension int, config *GraphSettings, matrix *sync.Pool, sk scoreKeeper, dynamicScore dynamicScoreKeeper, pool *sync.Pool) ([]cigar.Cigar, int64, int, int, []uint32) {
	s := pool.Get().(*dnaPool)
	defer pool.Put(s)

	s.Seq, s.Path = s.Seq[:0], s.Path[:0]
	s.Seq = getRightBases(n, extension, start, seq)
	s.Path = make([]uint32, len(currentPath))
	copy(s.Path, currentPath)
	if len(n.Seq)-start >= extension || len(n.Next) == 0 {
		sk.rightScore, sk.rightAlignment, _, sk.targetEnd, _, sk.queryEnd = RightLocal(s.Seq, read, config, matrix)
		sk.rightPath = s.Path
		
		return sk.rightAlignment, sk.rightScore, sk.targetEnd + start, sk.queryEnd, sk.rightPath
	} else {
		sk.rightScore = math.MinInt64
		for _, i := range n.Next {
			dynamicScore.route, s.currScore, s.targetEnd, s.queryEnd, s.Path = RightAlignTraversal(i.Dest, s.Seq, 0, s.Path, read, extension, config, matrix, sk, dynamicScore, pool)
			if s.currScore > sk.rightScore {
				sk.rightScore = s.currScore
				sk.rightAlignment = dynamicScore.route
				sk.targetEnd = s.targetEnd
				sk.queryEnd = s.queryEnd
				sk.rightPath = s.Path
			}
		}
		cigar.ReverseCigar(sk.rightAlignment)
		return sk.rightAlignment, sk.rightScore, sk.targetEnd + start, sk.queryEnd, sk.rightPath
	}
}

func ReversePath(alpha []uint32) {
	var i, off int
	for i = len(alpha)/2 - 1; i >= 0; i-- {
		off = len(alpha) - 1 - i
		alpha[i], alpha[off] = alpha[off], alpha[i]
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
