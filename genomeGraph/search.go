package genomeGraph

import (
	"bytes"
	"io"
	"log"
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

func getLeftBases(n *Node, extension int, start int, seq []dna.Base, ans []dna.Base) []dna.Base {
	currLen := len(seq)
	return append(append(ans, n.Seq[start-numbers.Min(currLen+start, extension)-currLen:start]...), seq...)
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

func LeftAlignTraversal(n *Node, seq []dna.Base, start int, currentPath []uint32, query []dna.Base, config *GraphSettings, sk scoreKeeper, memoryPool *sync.Pool) ([]cigar.Cigar, int64, int, int, []uint32) {
	alignData := memoryPool.Get().(*MemoryAllocation)
	defer memoryPool.Put(alignData)
	alignData.Alignment.Seq, alignData.Alignment.Path = alignData.Alignment.Seq[:0], alignData.Alignment.Path[:0]
	alignData.Alignment.Seq = getLeftBases(n, config.Extention, start, seq, alignData.Alignment.Seq)
	alignData.Alignment.Path = make([]uint32, len(currentPath))

	copy(alignData.Alignment.Path, currentPath)
	AddPath(alignData.Alignment.Path, n.Id)

	currLen := len(seq)
	if currLen+start >= config.Extention || len(n.Prev) == 0 {
		sk.leftScore, sk.leftAlignment, sk.targetStart, sk.queryStart = LeftDynamicAln(alignData.Alignment.Seq, query, config, alignData.Matrix)
		sk.targetStart = start - len(alignData.Alignment.Seq) - currLen + sk.targetStart
		sk.leftPath = alignData.Alignment.Path
		return sk.leftAlignment, sk.leftScore, sk.targetStart, sk.queryStart, sk.leftPath
	} else {
		//A very negative number
		sk.leftScore = math.MinInt64
		for _, i := range n.Prev {
			alignData.Matrix.route, alignData.Alignment.currScore, alignData.Alignment.targetStart, alignData.Alignment.queryStart, alignData.Alignment.Path = LeftAlignTraversal(i.Dest, alignData.Alignment.Seq, len(i.Dest.Seq), alignData.Alignment.Path, query, config, sk, memoryPool)
			if alignData.Alignment.currScore > sk.leftScore {
				sk.leftScore = alignData.Alignment.currScore
				sk.leftAlignment = alignData.Matrix.route
				sk.targetStart = start - len(alignData.Alignment.Seq) - currLen + alignData.Alignment.targetStart
				sk.queryStart = alignData.Alignment.queryStart
				sk.leftPath = alignData.Alignment.Path
			}
		}
		cigar.ReverseCigar(sk.leftAlignment)
		reversePath(sk.leftPath)
		return sk.leftAlignment, sk.leftScore, sk.targetStart, sk.queryStart, sk.leftPath
	}
}

func LeftDynamicAln(alpha []dna.Base, beta []dna.Base, config *GraphSettings, res *Matrix) (int64, []cigar.Cigar, int, int) {
	resetDynamicScore(res)
	columns, rows := len(alpha), len(beta)

	for res.i = 0; res.i <= columns; res.i++ {
		res.matrix[res.i][0] = 0
	}
	for res.j = 0; res.j <= rows; res.j++ {
		res.matrix[0][res.j] = 0
	}
	for res.i = 1; res.i <= columns; res.i++ {
		for res.j = 1; res.j <= rows; res.j++ {
			res.matrix[res.i][res.j], res.trace[res.i][res.j] = cigar.TripleMaxTrace(res.matrix[res.i-1][res.j-1]+config.ScoreMatrix[alpha[res.i-1]][beta[res.j-1]], res.matrix[res.i][res.j-1]+config.GapPenalty, res.matrix[res.i-1][res.j]+config.GapPenalty)
			res.matrix[res.i][res.j] = numbers.Max(res.matrix[res.i][res.j], 0)
		}
	}
	for res.i, res.j, res.index = columns, rows, 0; res.matrix[res.i][res.j] > 0; {
		if len(res.route) == 0 || res.route[len(res.route)-1].Op != res.trace[res.i][res.j] {
			res.route = append(res.route, cigar.Cigar{RunLength: 1, Op: res.trace[res.i][res.j]})
		} else {
			res.route[len(res.route)-1].RunLength++
		}
		switch res.trace[res.i][res.j] {
		case cigar.Match:
			res.i, res.j = res.i-1, res.j-1
		case cigar.Insertion:
			res.j -= 1
		case cigar.Deletion:
			res.i -= 1
		default:
			log.Fatalf("Error: unexpected traceback %c\n", res.trace[res.i][res.j])
		}
	}
	return res.matrix[columns][rows], res.route, res.i, res.j
}

func RightAlignTraversal(n *Node, seq []dna.Base, end int, currentPath []uint32, query []dna.Base, config *GraphSettings, sk scoreKeeper, memoryPool *sync.Pool) ([]cigar.Cigar, int64, int, int, []uint32) {
	alignData := memoryPool.Get().(*MemoryAllocation)
	defer memoryPool.Put(alignData)

	alignData.Alignment.Seq, alignData.Alignment.Path = alignData.Alignment.Seq[:0], alignData.Alignment.Path[:0]
	alignData.Alignment.Seq = getRightBases(n, config.Extention, end, seq, alignData.Alignment.Seq)
	alignData.Alignment.Path = make([]uint32, len(currentPath))
	copy(alignData.Alignment.Path, currentPath)

	if len(seq)+len(n.Seq)-end >= config.Extention || len(n.Next) == 0 {
		sk.rightScore, sk.rightAlignment, sk.targetEnd, sk.queryEnd = RightDynamicAln(alignData.Alignment.Seq, query, config, alignData.Matrix)
		sk.rightPath = alignData.Alignment.Path
		return sk.rightAlignment, sk.rightScore, sk.targetEnd + end, sk.queryEnd, sk.rightPath
	} else {
		sk.rightScore = math.MinInt64
		for _, i := range n.Next {
			alignData.Matrix.route, alignData.Alignment.currScore, alignData.Alignment.targetEnd, alignData.Alignment.queryEnd, alignData.Alignment.Path = RightAlignTraversal(i.Dest, alignData.Alignment.Seq, 0, alignData.Alignment.Path, query, config, sk, memoryPool)
			if alignData.Alignment.currScore > sk.rightScore {
				sk.rightScore = alignData.Alignment.currScore
				sk.rightAlignment = alignData.Matrix.route
				sk.targetEnd = alignData.Alignment.targetEnd
				sk.queryEnd = alignData.Alignment.queryEnd
				sk.rightPath = alignData.Alignment.Path
			}
		}
		cigar.ReverseCigar(sk.rightAlignment)
		return sk.rightAlignment, sk.rightScore, sk.targetEnd + end, sk.queryEnd, sk.rightPath
	}
}

func RightDynamicAln(alpha []dna.Base, beta []dna.Base, config *GraphSettings, res *Matrix) (int64, []cigar.Cigar, int, int) {
	columns, rows := len(alpha), len(beta)
	resetDynamicScore(res)
	var maxI int
	var maxJ int

	for res.i = 0; res.i <= columns; res.i++ {
		res.matrix[res.i][0] = config.GapPenalty * int64(res.i)
		res.trace[res.i][0] = cigar.Deletion // Initialize first column traces
	}
	for res.j = 0; res.j <= rows; res.j++ {
		res.matrix[0][res.j] = config.GapPenalty * int64(res.j)
		res.trace[0][res.j] = cigar.Insertion // Initialize first row traces
	}
	for res.i = 1; res.i <= columns; res.i++ {
		for res.j = 1; res.j <= rows; res.j++ {
			res.matrix[res.i][res.j], res.trace[res.i][res.j] = cigar.TripleMaxTrace(
				res.matrix[res.i-1][res.j-1]+config.ScoreMatrix[alpha[res.i-1]][beta[res.j-1]],
				res.matrix[res.i][res.j-1]+config.GapPenalty,
				res.matrix[res.i-1][res.j]+config.GapPenalty,
			)
			if res.matrix[res.i][res.j] > res.currMax {
				res.currMax = res.matrix[res.i][res.j]
				maxI = res.i
				maxJ = res.j
			}
		}
	}
	for res.i, res.j, res.index = maxI, maxJ, 0; res.i > 0 || res.j > 0; {
		if len(res.route) == 0 || res.route[len(res.route)-1].Op != res.trace[res.i][res.j] {
			res.route = append(res.route, cigar.Cigar{RunLength: 1, Op: res.trace[res.i][res.j]})
		} else {
			res.route[len(res.route)-1].RunLength++
		}
		switch res.trace[res.i][res.j] {
		case cigar.Match:
			res.i, res.j = res.i-1, res.j-1
		case cigar.Insertion:
			res.j -= 1
		case cigar.Deletion:
			res.i -= 1
		default:
			log.Fatalf("Error: unexpected traceback with %c\n", res.trace[res.i][res.j])
		}
	}
	return res.matrix[maxI][maxJ], res.route, maxI, maxJ
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
