package simpleGraph

import (
	"bytes"
	"github.com/vertgenlab/gonomics/cigar"
	"github.com/vertgenlab/gonomics/common"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/dnaTwoBit"
	"github.com/vertgenlab/gonomics/fastq"
	"github.com/vertgenlab/gonomics/fileio"
	"log"
	"strings"
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

type dynamicScoreKeeper struct {
	//scores [][]int64
	gapPen int64

	currMax int64
	route   []cigar.ByteCigar
	//maxI, maxJ int
	//minI, minJ int
}

func resetDynamicScore(sk dynamicScoreKeeper) {
	sk.route = sk.route[:0]
	sk.currMax = 0
}

func RightDynamicAln(alpha []dna.Base, beta []dna.Base, scores [][]int64, gapPen int64, m [][]int64, trace [][]byte, dynamicScore dynamicScoreKeeper) (int64, []cigar.ByteCigar, int, int, int, int) {
	//check if size of alpha is larger than m
	resetDynamicScore(dynamicScore)
	//var currMax int64
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
				m[i][j], trace[i][j] = cigar.ByteMatrixTrace(m[i-1][j-1]+scores[alpha[i-1]][beta[j-1]], m[i][j-1]+gapPen, m[i-1][j]+gapPen)
			}
			if m[i][j] > dynamicScore.currMax {
				dynamicScore.currMax = m[i][j]
				maxI = i
				maxJ = j
			}
		}
	}
	//var route []cigar.ByteCigar = make([]cigar.ByteCigar, 0, 1)
	//traceback starts in top corner
	curr := cigar.ByteCigar{}
	for i, j, routeIdx = maxI, maxJ, 0; i > 0 || j > 0; {
		//if route[routeIdx].RunLength == 0 {
		if len(dynamicScore.route) == 0 {
			curr = cigar.ByteCigar{RunLen: 1, Op: trace[i][j]}
			dynamicScore.route = append(dynamicScore.route, curr)
		} else if dynamicScore.route[routeIdx].Op == trace[i][j] {
			dynamicScore.route[routeIdx].RunLen += 1
		} else {
			curr = cigar.ByteCigar{RunLen: 1, Op: trace[i][j]}
			dynamicScore.route = append(dynamicScore.route, curr)
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
	cigar.ReverseBytesCigar(dynamicScore.route)
	return m[maxI][maxJ], dynamicScore.route, 0, maxI, 0, maxJ
}

func LeftDynamicAln(alpha []dna.Base, beta []dna.Base, scores [][]int64, gapPen int64, m [][]int64, trace [][]byte, dynamicScore dynamicScoreKeeper) (int64, []cigar.ByteCigar, int, int, int, int) {
	//check if size of alpha is larger than m
	resetDynamicScore(dynamicScore)
	var i, j, routeIdx int

	for i = 0; i < len(alpha)+1; i++ {
		m[i][0] = 0
	}
	for j = 0; j < len(beta)+1; j++ {
		m[0][j] = 0
	}
	for i = 1; i < len(alpha)+1; i++ {
		for j = 1; j < len(beta)+1; j++ {
			m[i][j], trace[i][j] = cigar.ByteMatrixTrace(m[i-1][j-1]+scores[alpha[i-1]][beta[j-1]], m[i][j-1]+gapPen, m[i-1][j]+gapPen)
			if m[i][j] < 0 {
				m[i][j] = 0
			}
		}
	}
	var minI, minJ = len(alpha), len(beta)
	//var route []cigar.ByteCigar = make([]cigar.ByteCigar, 0, 1)
	curr := cigar.ByteCigar{}
	//traceback starts in top corner
	for i, j, routeIdx = len(alpha), len(beta), 0; m[i][j] > 0; {
		//if route[routeIdx].RunLength == 0 {
		if len(dynamicScore.route) == 0 {
			curr = cigar.ByteCigar{RunLen: 1, Op: trace[i][j]}
			dynamicScore.route = append(dynamicScore.route, curr)
		} else if dynamicScore.route[routeIdx].Op == trace[i][j] {
			dynamicScore.route[routeIdx].RunLen += 1
		} else {
			curr = cigar.ByteCigar{RunLen: 1, Op: trace[i][j]}
			dynamicScore.route = append(dynamicScore.route, curr)
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
	cigar.ReverseBytesCigar(dynamicScore.route)
	return m[len(alpha)][len(beta)], dynamicScore.route, minI, len(alpha), minJ, len(beta)
}

func getRightBases(n *Node, ext int, start int, seq []dna.Base, ans []dna.Base) []dna.Base {
	var availableBases int = len(seq) + len(n.Seq) - start
	var targetLength int = common.Min(availableBases, ext)
	var basesToTake int = targetLength - len(seq)
	ans = make([]dna.Base, targetLength)
	copy(ans[0:len(seq)], seq)
	copy(ans[len(seq):targetLength], n.Seq[start:start+basesToTake])
	return ans
}

func RightAlignTraversal(n *Node, seq []dna.Base, start int, currentPath []uint32, ext int, read []dna.Base, m [][]int64, trace [][]byte, dnaPool *sync.Pool, dynamicScore dynamicScoreKeeper) ([]cigar.ByteCigar, int64, int, int, []uint32) {
	if len(seq) >= ext {
		log.Fatalf("Error: right traversal, the length=%d of DNA sequence in previous nodes should not be enough to satisfy the desired extenion=%d.\n", len(seq), ext)
	}

	s := dnaPool.Get().([]dna.Base)
	s = getRightBases(n, ext, start, seq, s)

	path := make([]uint32, len(currentPath))
	copy(path, currentPath)
	var bestTargetEnd, bestQueryEnd, targetEnd, queryEnd int
	var bestScore, score int64
	var bestAlignment, alignment []cigar.ByteCigar
	var bestPath []uint32

	if len(seq)+len(n.Seq)-start >= ext || len(n.Next) == 0 {
		score, alignment, _, targetEnd, _, queryEnd = RightDynamicAln(s, read, HumanChimpTwoScoreMatrix, -600, m, trace, dynamicScore)
		return alignment, score, targetEnd + start, queryEnd, path
	} else {
		bestScore = -1
		for _, i := range n.Next {
			alignment, score, targetEnd, queryEnd, path = RightAlignTraversal(i.Dest, s, 0, path, ext, read, m, trace, dnaPool, dynamicScore)
			if score > bestScore {
				bestScore = score
				bestAlignment = alignment
				bestTargetEnd = targetEnd
				bestQueryEnd = queryEnd
				bestPath = path
			}
		}

	}
	s = s[:0]
	dnaPool.Put(s)
	return bestAlignment, bestScore, bestTargetEnd + start, bestQueryEnd, bestPath
}

func getLeftTargetBases(n *Node, ext int, refEnd int, seq []dna.Base, s []dna.Base) []dna.Base {
	var availableBases int = len(seq) + refEnd
	var targetLength int = common.Min(availableBases, ext)
	var basesToTake int = targetLength - len(seq)
	s = make([]dna.Base, targetLength)
	copy(s[0:basesToTake], n.Seq[refEnd-basesToTake:refEnd])
	copy(s[basesToTake:targetLength], seq)
	return s
}

func LeftAlignTraversal(n *Node, seq []dna.Base, refEnd int, currentPath []uint32, ext int, read []dna.Base, m [][]int64, trace [][]byte, dnaPool *sync.Pool, dynamicScore dynamicScoreKeeper) ([]cigar.ByteCigar, int64, int, int, []uint32) {
	if len(seq) >= ext {
		log.Fatalf("Error: left traversal, the length=%d of DNA sequence in previous nodes should not be enough to satisfy the desired extenion=%d.\n", len(seq), ext)
	}
	s := dnaPool.Get().([]dna.Base)
	s = getLeftTargetBases(n, ext, refEnd, seq, s)

	path := make([]uint32, len(currentPath))
	copy(path, currentPath)
	AddPath(path, n.Id)
	var bestQueryStart, queryStart, refStart, bestRefStart int
	var bestScore, score int64
	var bestAlignment, alignment []cigar.ByteCigar
	var bestPath []uint32

	if len(seq)+refEnd >= ext || len(n.Prev) == 0 {
		score, alignment, refStart, _, queryStart, _ = LeftDynamicAln(s, read, HumanChimpTwoScoreMatrix, -600, m, trace, dynamicScore)

		refEnd = refEnd - len(s) - len(seq) + refStart
		return alignment, score, refEnd, queryStart, path
	} else {
		bestScore = -1
		for _, i := range n.Prev {
			alignment, score, refStart, queryStart, path = LeftAlignTraversal(i.Dest, s, len(i.Dest.Seq), path, ext, read, m, trace, dnaPool, dynamicScore)
			if score > bestScore {
				bestScore = score
				bestAlignment = alignment
				bestRefStart = refEnd - len(s) - len(seq) + refStart
				bestQueryStart = queryStart
				bestPath = path
			}
		}
	}
	s = s[:0]
	dnaPool.Put(s)
	return bestAlignment, bestScore, bestRefStart, bestQueryStart, bestPath
}

func ReversePath(alpha []uint32) {
	var i, off int
	for i, off = len(alpha)/2-1, len(alpha)-1-i; i >= 0; i-- {
		alpha[i], alpha[off] = alpha[off], alpha[i]
	}
}

func readFastqGsw(fileOne string, fileTwo string, answer chan<- *fastq.PairedEndBig) {
	readOne, readTwo := fileio.NewSimpleReader(fileOne), fileio.NewSimpleReader(fileTwo)
	for fq, done := fqPair(readOne, readTwo); !done; fq, done = fqPair(readOne, readTwo) {
		answer <- fq
	}
	close(answer)
}

func fqPair(reader1 *fileio.SimpleReader, reader2 *fileio.SimpleReader) (*fastq.PairedEndBig, bool) {
	fqOne, done1 := nextFq(reader1)
	fqTwo, done2 := nextFq(reader2)
	if (!done1 && done2) || (done1 && !done2) {
		log.Fatalf("Error: fastq files do not end at the same time...\n")
	} else if done1 || done2 {
		return nil, true
	}
	curr := NewBigFastqPair()
	curr.Fwd, curr.Rev = fqOne, fqTwo
	//curr.Fwd.Name, curr.Rev.Name = strings.Split(fqOne.Name, " ")[0], strings.Split(fqTwo.Name, " ")[0]
	return curr, false
}

func nextFq(reader *fileio.SimpleReader) (*fastq.FastqBig, bool) {
	name, done := fileio.ReadLine(reader)
	if done {
		return nil, true
	}
	seq, sDone := fileio.ReadLine(reader)
	plus, pDone := fileio.ReadLine(reader)
	qual, qDone := fileio.ReadLine(reader)

	if sDone || pDone || qDone {
		log.Fatalf("Error: There is an empty line in this fastq record\n")
	}

	answer := fastq.FastqBig{}
	//data := simplePool.Get().(*bytes.Buffer)
	//data.Write(name)
	answer.Name = strings.Split(string(name[1:]), " ")[0]
	//data.Reset()
	//data.Write(seq)
	//set up sequence and reverse comp
	answer.Seq = ByteSliceToDnaBases(seq)
	answer.SeqRc = make([]dna.Base, len(answer.Seq))
	copy(answer.SeqRc, answer.Seq)
	dna.ReverseComplement(answer.SeqRc)

	//performs two bit conversion
	answer.Rainbow = dnaTwoBit.NewTwoBitRainbow(answer.Seq)
	answer.RainbowRc = dnaTwoBit.NewTwoBitRainbow(answer.SeqRc)

	//data.Reset()

	if string(plus) != "+" {
		log.Fatalf("Error: This line should be a + (plus) sign \n")
	}
	//data.Write(qual)
	answer.Qual = fastq.ToQualUint8(bytes.Runes(qual))

	//data.Reset()
	//simplePool.Put(data)
	return &answer, false
}

func NewBigFastqPair() *fastq.PairedEndBig {
	return &fastq.PairedEndBig{
		Fwd: new(fastq.FastqBig),
		Rev: new(fastq.FastqBig),
	}
}

func ByteToBase(b byte) dna.Base {
	switch b {
	case 'A':
		return dna.A
	case 'C':
		return dna.C
	case 'G':
		return dna.G
	case 'T':
		return dna.T
	case 'N':
		return dna.N
	case 'a':
		return dna.A
	case 'c':
		return dna.C
	case 'g':
		return dna.G
	case 't':
		return dna.T
	case 'n':
		return dna.N
	case '-':
		return dna.Gap
	//VCF uses star to denote a deleted allele
	case '*':
		return dna.Gap
	case '.':
		return dna.Dot
	default:
		log.Fatalf("Error: unexpected character in dna %c\n", b)
		return dna.N
	}
}

func ByteSliceToDnaBases(b []byte) []dna.Base {
	var answer []dna.Base = make([]dna.Base, len(b))
	for i, byteValue := range b {
		answer[i] = ByteToBase(byteValue)
	}
	return answer
}
