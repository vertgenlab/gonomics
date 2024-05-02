package genomeGraph

import (
	"log"
	"math"
	"sync"

	"github.com/vertgenlab/gonomics/cigar"
	"github.com/vertgenlab/gonomics/dna/dnaTwoBit"

	"github.com/vertgenlab/gonomics/numbers"
)

func LeftGetTwoBit(n *Node, extension int, position int, seq *dnaTwoBit.TwoBit, ans *dnaTwoBit.TwoBit) *dnaTwoBit.TwoBit {
	basesToTake := numbers.Min(seq.Len+position, extension) - seq.Len
	if basesToTake > 0 {
		frag := dnaTwoBit.GetFrag(n.SeqTwoBit, position-basesToTake, position)
		dnaTwoBit.Cat(ans, frag)
	}
	dnaTwoBit.Cat(ans, seq)
	return ans
}

func TwoBitLocalLeftAlign(alpha *dnaTwoBit.TwoBit, beta *dnaTwoBit.TwoBit, settings *GraphSettings, pool *sync.Pool) (int64, []cigar.ByteCigar, int, int) {
	rows, columns := alpha.Len+1, beta.Len+1
	matrix := pool.Get().(*MatrixMemoryPool)
	matrix.Route = matrix.Route[:0]
	matrix.CurrMax = 0

	if cap(matrix.matrix) < rows || cap(matrix.matrix[0]) < columns {
		matrix.matrix = make([][]int64, rows)
		matrix.trace = make([][]byte, rows)
		for idx := range matrix.matrix {
			matrix.matrix[idx] = make([]int64, columns)
			matrix.trace[idx] = make([]byte, columns)
		}
	}
	for matrix.i = 0; matrix.i < rows; matrix.i++ {
		matrix.matrix[matrix.i][0] = 0
	}

	for matrix.j = 0; matrix.j < columns; matrix.j++ {
		matrix.matrix[0][matrix.j] = 0
	}

	for matrix.i = 1; matrix.i < rows; matrix.i++ {
		for matrix.j = 1; matrix.j < columns; matrix.j++ {
			matrix.matrix[matrix.i][matrix.j], matrix.trace[matrix.i][matrix.j] = cigar.TraceMatrixExtension(matrix.matrix[matrix.i-1][matrix.j-1], matrix.matrix[matrix.i-1][matrix.j-1]+settings.ScoreMatrix[uint(dnaTwoBit.GetBase(alpha, uint(matrix.i-1)))][uint(dnaTwoBit.GetBase(beta, uint(matrix.j-1)))], matrix.matrix[matrix.i][matrix.j-1]+settings.GapPenalty, matrix.matrix[matrix.i-1][matrix.j]+settings.GapPenalty)
			if matrix.matrix[matrix.i][matrix.j] < 0 {
				matrix.matrix[matrix.i][matrix.j] = 0
			}
		}
	}

	//traceback starts in top corner
	var minI, minJ int
	for matrix.i, matrix.j, matrix.routeIdx = rows-1, columns-1, 0; matrix.matrix[matrix.i][matrix.j] > 0; {
		//if route[matrix.routeIdx].RunLength == 0 {
		if len(matrix.Route) == 0 {
			curr := cigar.ByteCigar{RunLen: 1, Op: matrix.trace[matrix.i][matrix.j]}
			matrix.Route = append(matrix.Route, curr)
		} else if matrix.Route[matrix.routeIdx].Op == matrix.trace[matrix.i][matrix.j] {
			matrix.Route[matrix.routeIdx].RunLen += 1
		} else {
			curr := cigar.ByteCigar{RunLen: 1, Op: matrix.trace[matrix.i][matrix.j]}
			matrix.Route = append(matrix.Route, curr)
			matrix.routeIdx++
		}
		switch matrix.trace[matrix.i][matrix.j] {
		case '=':
			matrix.i, matrix.j = matrix.i-1, matrix.j-1
		case 'X':
			matrix.i, matrix.j = matrix.i-1, matrix.j-1
		case 'I':
			matrix.j -= 1
		case 'D':
			matrix.i -= 1
		default:
			log.Fatalf("Error: unexpected traceback %c\n", matrix.trace[matrix.i][matrix.j])
		}
		minI = matrix.i
		minJ = matrix.j
	}
	cigar.ReverseBytesCigar(matrix.Route)
	return matrix.matrix[rows-1][columns-1], matrix.Route, minI, minJ
}

func LeftTwoBitDfs(n *Node, seq *dnaTwoBit.TwoBit, position int, currentPath []uint32, read *dnaTwoBit.TwoBit, settings *GraphSettings, sk scoreKeeper, memory *sync.Pool) ([]cigar.ByteCigar, int64, int, int, []uint32) {

	s := memory.Get().(*MatrixMemoryPool)
	defer memory.Put(s)

	//currSeq, s.Path = currSeq, s.Path[:0]
	var currSeq *dnaTwoBit.TwoBit
	currSeq = LeftGetTwoBit(n, settings.Extension, position, seq, currSeq)
	s.Path = append(s.Path[:0], currentPath...)
	AddPath(s.Path, n.Id)

	if seq.Len+position >= settings.Extension || len(n.Prev) == 0 {

		sk.leftScore, sk.leftAlignment, sk.targetStart, sk.queryStart = TwoBitLocalLeftAlign(currSeq, read, settings, memory)
		sk.targetStart = position - currSeq.Len - seq.Len + sk.targetStart
		sk.leftPath = s.Path
		return sk.leftAlignment, sk.leftScore, sk.targetStart, sk.queryStart, sk.leftPath
	}

	sk.leftScore = math.MinInt64
	for _, i := range n.Prev {
		s.Route, s.CurrScore, s.TargetStart, s.QueryStart, s.Path = LeftTwoBitDfs(i.Dest, currSeq, len(i.Dest.Seq), s.Path, read, settings, sk, memory)
		if s.CurrScore > sk.leftScore {
			sk.leftScore = s.CurrScore
			sk.leftAlignment = s.Route
			sk.targetStart = s.TargetStart
			sk.queryStart = s.QueryStart
			sk.leftPath = s.Path
		}
	}

	ReversePath(sk.leftPath)
	return sk.leftAlignment, sk.leftScore, sk.targetStart, sk.queryStart, sk.leftPath
}

func RightGetTwoBit(n *Node, extension int, position int, seq *dnaTwoBit.TwoBit, ans *dnaTwoBit.TwoBit) *dnaTwoBit.TwoBit {
	basesToTake := numbers.Min(n.SeqTwoBit.Len-position, extension)
	dnaTwoBit.Cat(ans, seq)
	if basesToTake > 0 {
		frag := dnaTwoBit.GetFrag(n.SeqTwoBit, position, position+basesToTake)
		if frag != nil {
			dnaTwoBit.Cat(ans, frag)
		}
	}
	return ans
}

func TwoBitLocalRightAlign(alpha *dnaTwoBit.TwoBit, beta *dnaTwoBit.TwoBit, settings *GraphSettings, pool *sync.Pool) (int64, []cigar.ByteCigar, int, int) {
	rows, columns := alpha.Len+1, beta.Len+1
	matrix := pool.Get().(*MatrixMemoryPool)
	defer pool.Put(matrix)

	var maxI int
	var maxJ int
	var currMax int64

	if cap(matrix.matrix) < rows || cap(matrix.matrix[0]) < columns {
		matrix.matrix = make([][]int64, rows)
		matrix.trace = make([][]byte, rows)
		for idx := range matrix.matrix {
			matrix.matrix[idx] = make([]int64, columns)
			matrix.trace[idx] = make([]byte, columns)
		}
	}

	for matrix.i = 0; matrix.i < rows; matrix.i++ {
		for matrix.j = 0; matrix.j < columns; matrix.j++ {
			if matrix.i == 0 && matrix.j == 0 {
				matrix.matrix[matrix.i][matrix.j] = 0
			} else if matrix.i == 0 {
				matrix.matrix[matrix.i][matrix.j] = matrix.matrix[matrix.i][matrix.j-1] + settings.GapPenalty
				matrix.trace[matrix.i][matrix.j] = cigar.Insertion
			} else if matrix.j == 0 {
				matrix.matrix[matrix.i][matrix.j] = matrix.matrix[matrix.i-1][matrix.j] + settings.GapPenalty
				matrix.trace[matrix.i][matrix.j] = cigar.Deletion
			} else {
				matrix.matrix[matrix.i][matrix.j], matrix.trace[matrix.i][matrix.j] = cigar.TraceMatrixExtension(matrix.matrix[matrix.i-1][matrix.j-1], matrix.matrix[matrix.i-1][matrix.j-1]+settings.ScoreMatrix[uint(dnaTwoBit.GetBase(alpha, uint(matrix.i-1)))][uint(dnaTwoBit.GetBase(beta, uint(matrix.j-1)))], matrix.matrix[matrix.i][matrix.j-1]+settings.GapPenalty, matrix.matrix[matrix.i-1][matrix.j]+settings.GapPenalty)
			}
			if matrix.matrix[matrix.i][matrix.j] > currMax {
				currMax = matrix.matrix[matrix.i][matrix.j]
				maxI = matrix.i
				maxJ = matrix.j
			}
		}
	}
	//traceback starts in top corner
	for matrix.i, matrix.j, matrix.routeIdx = maxI, maxJ, 0; matrix.i > 0 || matrix.j > 0; {
		//if route[matrix.routeIdx].RunLength == 0 {
		if len(matrix.Route) == 0 {
			curr := cigar.ByteCigar{RunLen: 1, Op: matrix.trace[matrix.i][matrix.j]}
			matrix.Route = append(matrix.Route, curr)
		} else if matrix.Route[matrix.routeIdx].Op == matrix.trace[matrix.i][matrix.j] {
			matrix.Route[matrix.routeIdx].RunLen += 1
		} else {
			curr := cigar.ByteCigar{RunLen: 1, Op: matrix.trace[matrix.i][matrix.j]}
			matrix.Route = append(matrix.Route, curr)
			matrix.routeIdx++
		}
		switch matrix.trace[matrix.i][matrix.j] {
		case cigar.Equal:
			matrix.i, matrix.j = matrix.i-1, matrix.j-1
		case cigar.Mismatch:
			matrix.i, matrix.j = matrix.i-1, matrix.j-1
		case cigar.Insertion:
			matrix.j -= 1
		case cigar.Deletion:
			matrix.i -= 1
		default:
			log.Fatalf("Error: unexpected traceback with %c\n", matrix.trace[matrix.i][matrix.j])
		}
	}
	cigar.ReverseBytesCigar(matrix.Route)
	return matrix.matrix[maxI][maxJ], matrix.Route, maxI, maxJ
}

func RightTwoBitDfs(n *Node, seq *dnaTwoBit.TwoBit, start int, currentPath []uint32, read *dnaTwoBit.TwoBit, settings *GraphSettings, sk *scoreKeeper, memory *sync.Pool) ([]cigar.ByteCigar, int64, int, int, []uint32) {
	s := memory.Get().(*MatrixMemoryPool)
	defer memory.Put(s)
	var currSeq *dnaTwoBit.TwoBit

	currSeq, s.Path = RightGetTwoBit(n, settings.Extension, start, seq, currSeq), append([]uint32(nil), currentPath...) // Reuse memory more efficiently
	if seq.Len+n.SeqTwoBit.Len-start >= settings.Extension || len(n.Next) == 0 {
		sk.rightScore, sk.rightAlignment, sk.targetEnd, sk.queryEnd = TwoBitLocalRightAlign(currSeq, read, settings, memory)
		sk.rightPath = s.Path
		return sk.rightAlignment, sk.rightScore, sk.targetEnd + start, sk.queryEnd, sk.rightPath
	}

	sk.rightScore = math.MinInt64
	for _, i := range n.Next {
		alignment, score, targetEnd, queryEnd, path := LeftTwoBitDfs(i.Dest, currSeq, 0, s.Path, read, settings, *sk, memory)
		if score > sk.rightScore {
			sk.rightScore, sk.rightAlignment, sk.targetEnd, sk.queryEnd, sk.rightPath = score, alignment, targetEnd, queryEnd, path
		}
	}
	cigar.ReverseBytesCigar(sk.rightAlignment)
	return sk.rightAlignment, sk.rightScore, sk.targetEnd + start, sk.queryEnd, sk.rightPath
}
