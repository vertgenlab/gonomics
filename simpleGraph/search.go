package simpleGraph

import (
	"github.com/vertgenlab/gonomics/cigar"
	"github.com/vertgenlab/gonomics/common"
	"github.com/vertgenlab/gonomics/dna"
	"log"
	"sync"
	"math"
)

func NewMemSeedPool() sync.Pool {
	return sync.Pool{
		New: func() interface{} {
			pool := memoryPool{
				Hits:   make([]*SeedDev, 0, 10000),
				Worker: make([]*SeedDev, 0, 10000),
			}
			return &pool
		},
	}
}

type memoryPool struct {
	Hits   []*SeedDev
	Worker []*SeedDev
}

/*
func DevGraphSmithWaterman(gg *SimpleGraph, read fastq.FastqBig, seedHash map[uint64][]uint64, seedLen int, stepSize int, matrix *MatrixAln, scoreMatrix [][]int64, seedPool *sync.Pool, dnaPool *sync.Pool, sk scoreKeeper, dynamicScore dynamicScoreKeeper) giraf.Giraf {
	var currBest giraf.Giraf = giraf.Giraf{
		QName:     read.Name,
		QStart:    0,
		QEnd:      0,
		PosStrand: true,
		Path:      giraf.Path{},
		Cigar:     nil,
		AlnScore:  0,
		MapQ:      255,
		Seq:       read.Seq,
		Qual:      read.Qual,
		Notes:     []giraf.Note{giraf.Note{Tag: []byte("XO")[:2], Type: 'Z', Value: "~"}},
	}
	resetScoreKeeper(sk)
	sk.perfectScore = perfectMatchBig(&read, scoreMatrix)
	sk.extension = int(sk.perfectScore/600) + len(read.Seq)
	seeds := seedPool.Get().(*memoryPool)
	seeds.Hits = seeds.Hits[:0]
	seeds.Worker = seeds.Worker[:0]
	//tempSeeds := extendSeeds.Get().([]*SeedDev)
	//var tmpSeeds []*SeedDev = extendPool.Get().([]*SeedDev)
	seeds.Hits = seedMapMemPool(seedHash, gg.Nodes, &read, seedLen, sk.perfectScore, scoreMatrix, seeds.Hits, seeds.Worker)

	//tempSeeds = tempSeeds[:0]
	//extendSeeds.Put(tempSeeds)
	SortSeedDevByLen(seeds.Hits)
	var tailSeed *SeedDev
	//var seedScore int64
	var currSeq []dna.Base = make([]dna.Base, len(read.Seq))
	var currSeed *SeedDev

	//log.Printf("Number of hits: %d\n", len(seeds.Hits))
	for i := 0; i < len(seeds.Hits) && seedCouldBeBetter(int64(seeds.Hits[i].TotalLength), int64(currBest.AlnScore), sk.perfectScore, int64(len(read.Seq)), 100, 90, -196, -296); i++ {
		currSeed = seeds.Hits[i]
		tailSeed = getLastPart(currSeed)
		if currSeed.PosStrand {
			currSeq = read.Seq
		} else {
			currSeq = read.SeqRc
		}
		sk.seedScore = scoreSeedSeq(currSeq, currSeed.QueryStart, tailSeed.QueryStart+tailSeed.Length, scoreMatrix)
		if int(currSeed.TotalLength) == len(currSeq) {
			sk.minTarget = int(currSeed.TargetStart)
			sk.maxTarget = int(tailSeed.TargetStart + tailSeed.Length)
			sk.minQuery = int(currSeed.QueryStart)
			//sk.maxQuery = int(currSeed.TotalLength - 1)
			sk.currScore = sk.seedScore
		} else {

			sk.leftAlignment, sk.leftScore, sk.minTarget, sk.minQuery, sk.leftPath = LeftAlignTraversal(gg.Nodes[currSeed.TargetId], sk.leftSeq, int(currSeed.TargetStart), sk.leftPath, sk.extension-int(currSeed.TotalLength), currSeq[:currSeed.QueryStart], scoreMatrix, matrix, dynamicScore, dnaPool)
			sk.rightAlignment, sk.rightScore, sk.maxTarget, sk.maxQuery, sk.rightPath = RightAlignTraversal(gg.Nodes[tailSeed.TargetId], sk.rightSeq, int(tailSeed.TargetStart+tailSeed.Length), sk.rightPath, sk.extension-int(currSeed.TotalLength), currSeq[tailSeed.QueryStart+tailSeed.Length:], scoreMatrix, matrix, dynamicScore, dnaPool)
			sk.currScore = sk.leftScore + sk.seedScore + sk.rightScore

		}
		if sk.currScore > int64(currBest.AlnScore) {
			currBest.QStart = sk.minQuery
			currBest.QEnd = int(currSeed.QueryStart) + sk.minQuery + sk.maxQuery + int(currSeed.TotalLength) - 1
			currBest.PosStrand = currSeed.PosStrand
			ReversePath(sk.leftPath)
			currBest.Path = setPath(currBest.Path, sk.minTarget, CatPaths(CatPaths(sk.leftPath, getSeedPath(currSeed)), sk.rightPath), sk.maxTarget)
			currBest.Cigar = SoftClipBases(sk.minQuery, len(currSeq), cigar.CatByteCigar(cigar.AddCigarByte(sk.leftAlignment, cigar.ByteCigar{RunLen: uint16(sumLen(currSeed)), Op: 'M'}), sk.rightAlignment))
			currBest.AlnScore = int(sk.currScore)
			currBest.Seq = currSeq
			if &gg.Nodes[currBest.Path.Nodes[0]].Info != nil {
				currBest.Notes[0].Value = fmt.Sprintf("%s=%d", gg.Nodes[currBest.Path.Nodes[0]].Name, gg.Nodes[currBest.Path.Nodes[0]].Info.Start)
				currBest.Notes = append(currBest.Notes, infoToNotes(gg.Nodes, currBest.Path.Nodes))
			} else {
				currBest.Notes[0].Value = fmt.Sprintf("%s=%d", gg.Nodes[currBest.Path.Nodes[0]].Name, 1)
			}
		}
	}
	seedPool.Put(seeds)
	if !currBest.PosStrand {
		fastq.ReverseQualUint8Record(currBest.Qual)
	}
	return currBest
}

func DevWrapPairGiraf(gg *SimpleGraph, fq fastq.PairedEndBig, seedHash map[uint64][]uint64, seedLen int, stepSize int, matrix *MatrixAln, scoreMatrix [][]int64, seedPool *sync.Pool, dnaPool *sync.Pool, sk scoreKeeper, dynamicScore dynamicScoreKeeper) giraf.GirafPair {
	var mappedPair giraf.GirafPair = giraf.GirafPair{
		ReadOne: DevGraphSmithWaterman(gg, fq.Fwd, seedHash, seedLen, stepSize, matrix, scoreMatrix, seedPool, dnaPool, sk, dynamicScore),
		ReadTwo: DevGraphSmithWaterman(gg, fq.Rev, seedHash, seedLen, stepSize, matrix, scoreMatrix, seedPool, dnaPool, sk, dynamicScore),
	}
	//setGirafFlags(&mappedPair)
	return mappedPair
}*/

type MatrixAln struct {
	m     [][]int64
	trace [][]byte
}

type dnaPool struct {
	Seq                  []dna.Base
	Path                 []uint32
	queryStart, refStart int
	targetEnd, queryEnd  int
}

func NewDnaPool() sync.Pool {
	return sync.Pool{
		New: func() interface{} {
			dnaSeq := dnaPool{
				Seq:        make([]dna.Base, 0, 150),
				Path:       make([]uint32, 0, 10),
				queryStart: 0,
				refStart:   0,
				targetEnd:  0,
				queryEnd:   0,
			}
			return &dnaSeq
		},
	}
}

type dynamicScoreKeeper struct {
	gapPen  int64
	currMax int64
	route   []cigar.ByteCigar
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

type scoreKeeper struct {
	minTarget      int
	maxTarget      int
	minQuery       int
	maxQuery       int
	extension      int
	currScore      int64
	seedScore      int64
	perfectScore   int64
	leftScore      int64
	rightScore     int64
	leftPath       []uint32
	rightPath      []uint32
	leftSeq        []dna.Base
	rightSeq       []dna.Base
	currSeq        []dna.Base
	leftAlignment  []cigar.ByteCigar
	rightAlignment []cigar.ByteCigar
}

func resetScoreKeeper(sk scoreKeeper) {
	sk.minTarget, sk.maxTarget = 0, 0
	sk.minQuery, sk.maxQuery = 0, 0
	sk.currScore, sk.seedScore = 0, 0
	sk.perfectScore = 0
	sk.leftAlignment, sk.rightAlignment = sk.leftAlignment[:0], sk.rightAlignment[:0]
	sk.leftPath, sk.rightPath = sk.leftPath[:0], sk.rightPath[:0]
	sk.leftSeq, sk.rightSeq = sk.leftSeq[:0], sk.rightSeq[:0]
	sk.leftScore, sk.rightScore = 0, 0
	sk.currSeq = sk.currSeq[:0]
}

func getRightBases(n *Node, ext int, start int, seq []dna.Base, ans []dna.Base) []dna.Base {
	var availableBases int = len(seq) + len(n.Seq) - start
	var targetLength int = common.Min(availableBases, ext)
	var basesToTake int = targetLength - len(seq)
	return append(append(ans, seq...), n.Seq[start:start+basesToTake]...)
}
/*
func rightBasesFromTwoBit(n *Node, ext int, start int, seq []dna.Base, ans []dna.Base) []dna.Base {
	var availableBases int = len(seq) + len(n.Seq) - start
	var targetLength int = common.Min(availableBases, ext)
	var basesToTake int = targetLength - len(seq)

	return append(append(ans, seq...), dnaTwoBit.GetFrag(n.SeqTwoBit, start, start+basesToTake)...)
}*/

func RightAlignTraversal(n *Node, seq []dna.Base, start int, currentPath []uint32, ext int, read []dna.Base, scoreMatrix [][]int64, matrix *MatrixAln, sk scoreKeeper, dynamicScore dynamicScoreKeeper, pool *sync.Pool) ([]cigar.ByteCigar, int64, int, int, []uint32) {
	if len(seq) >= ext {
		log.Fatalf("Error: right traversal, the length=%d of DNA sequence in previous nodes should not be enough to satisfy the desired extenion=%d.\n", len(seq), ext)
	}
	s := pool.Get().(*dnaPool)
	s.Seq = s.Seq[:0]
	s.Seq = getRightBases(n, ext, start, seq, s.Seq)
	s.Path = s.Path[:0]
	s.Path = append(s.Path, currentPath...)
	if len(seq)+len(n.Seq)-start >= ext || len(n.Next) == 0 {
		sk.rightScore, sk.rightAlignment, _, sk.maxTarget, _, sk.maxQuery = RightDynamicAln(s.Seq, read, scoreMatrix, matrix, -600, dynamicScore)
		sk.rightPath = s.Path
		pool.Put(s)
		return sk.rightAlignment, sk.rightScore, sk.maxTarget + start, sk.maxQuery, sk.rightPath
	} else {
		sk.rightScore = math.MinInt64
		var score int64
		for _, i := range n.Next {
			dynamicScore.route, score, s.targetEnd, s.queryEnd, s.Path = RightAlignTraversal(i.Dest, s.Seq, 0, s.Path, ext, read, scoreMatrix, matrix, sk, dynamicScore, pool)
			if score > sk.rightScore {
				sk.rightScore = score
				sk.rightAlignment = dynamicScore.route
				sk.maxTarget = s.targetEnd
				sk.maxQuery = s.queryEnd
				sk.rightPath = s.Path
			}
		}
		pool.Put(s)
		cigar.ReverseBytesCigar(sk.rightAlignment)
		return sk.rightAlignment, sk.rightScore, sk.maxTarget + start, sk.maxQuery, sk.rightPath
	}
}

func RightDynamicAln(alpha []dna.Base, beta []dna.Base, scores [][]int64, matrix *MatrixAln, gapPen int64, dynamicScore dynamicScoreKeeper) (int64, []cigar.ByteCigar, int, int, int, int) {
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
				matrix.m[i][j] = 0
			} else if i == 0 {
				matrix.m[i][j] = matrix.m[i][j-1] + gapPen
				matrix.trace[i][j] = 'I'
			} else if j == 0 {
				matrix.m[i][j] = matrix.m[i-1][j] + gapPen
				matrix.trace[i][j] = 'D'
			} else {
				matrix.m[i][j], matrix.trace[i][j] = cigar.ByteMatrixTrace(matrix.m[i-1][j-1]+scores[alpha[i-1]][beta[j-1]], matrix.m[i][j-1]+gapPen, matrix.m[i-1][j]+gapPen)
			}
			if matrix.m[i][j] > dynamicScore.currMax {
				dynamicScore.currMax = matrix.m[i][j]
				maxI = i
				maxJ = j
			}
		}

	}

	//var route []cigar.ByteCigar = make([]cigar.ByteCigar, 0, 1)
	//traceback starts in top corner
	//dynamicScore.curr := cigar.ByteCigar{}
	for i, j, routeIdx = maxI, maxJ, 0; i > 0 || j > 0; {
		//if route[routeIdx].RunLength == 0 {
		if len(dynamicScore.route) == 0 {
			dynamicScore.route = append(dynamicScore.route, cigar.ByteCigar{RunLen: 1, Op: matrix.trace[i][j]})
		} else if dynamicScore.route[routeIdx].Op == matrix.trace[i][j] {
			dynamicScore.route[routeIdx].RunLen += 1
		} else {
			//dynamicScore.curr =
			dynamicScore.route = append(dynamicScore.route, cigar.ByteCigar{RunLen: 1, Op: matrix.trace[i][j]})
			routeIdx++
		}
		switch matrix.trace[i][j] {
		case 'M':
			i, j = i-1, j-1
		case 'I':
			j -= 1
		case 'D':
			i -= 1
		default:
			log.Fatalf("Error: unexpected traceback with %c\n", matrix.trace[i][j])
		}
	}
	//cigar.ReverseBytesCigar(dynamicScore.route)
	return matrix.m[maxI][maxJ], dynamicScore.route, 0, maxI, 0, maxJ
}

func LeftDynamicAln(alpha []dna.Base, beta []dna.Base, scores [][]int64, matrix *MatrixAln, gapPen int64, dynamicScore dynamicScoreKeeper) (int64, []cigar.ByteCigar, int, int, int, int) {
	//check if size of alpha is larger than m
	resetDynamicScore(dynamicScore)
	var i, j, routeIdx int

	for i = 0; i < len(alpha)+1; i++ {
		matrix.m[i][0] = 0
	}
	for j = 0; j < len(beta)+1; j++ {
		matrix.m[0][j] = 0
	}
	for i = 1; i < len(alpha)+1; i++ {
		for j = 1; j < len(beta)+1; j++ {
			matrix.m[i][j], matrix.trace[i][j] = cigar.ByteMatrixTrace(matrix.m[i-1][j-1]+scores[alpha[i-1]][beta[j-1]], matrix.m[i][j-1]+gapPen, matrix.m[i-1][j]+gapPen)
			if matrix.m[i][j] < 0 {
				matrix.m[i][j] = 0
			}
		}
	}
	var minI, minJ = len(alpha), len(beta)
	//var route []cigar.ByteCigar = make([]cigar.ByteCigar, 0, 1)
	//curr := cigar.ByteCigar{}
	//traceback starts in top corner
	for i, j, routeIdx = len(alpha), len(beta), 0; matrix.m[i][j] > 0; {
		//if route[routeIdx].RunLength == 0 {
		if len(dynamicScore.route) == 0 {
			dynamicScore.route = append(dynamicScore.route, cigar.ByteCigar{RunLen: 1, Op: matrix.trace[i][j]})
		} else if dynamicScore.route[routeIdx].Op == matrix.trace[i][j] {
			dynamicScore.route[routeIdx].RunLen += 1
		} else {
			dynamicScore.route = append(dynamicScore.route, cigar.ByteCigar{RunLen: 1, Op: matrix.trace[i][j]})
			routeIdx++
		}
		switch matrix.trace[i][j] {
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
	//cigar.ReverseBytesCigar(dynamicScore.route)
	return matrix.m[len(alpha)][len(beta)], dynamicScore.route, minI, len(alpha), minJ, len(beta)
}

func getLeftTargetBases(n *Node, ext int, refEnd int, seq []dna.Base, ans []dna.Base) []dna.Base {
	var availableBases int = len(seq) + refEnd
	var targetLength int = common.Min(availableBases, ext)
	var basesToTake int = targetLength - len(seq)
	return append(append(ans, n.Seq[refEnd-basesToTake:refEnd]...), seq...)
}
/*
func leftBasesFromTwoBit(n *Node, ext int, refEnd int, seq []dna.Base, ans []dna.Base) []dna.Base {
	var availableBases int = len(seq) + refEnd
	var targetLength int = common.Min(availableBases, ext)
	var basesToTake int = targetLength - len(seq)
	return append(append(ans, seq...), dnaTwoBit.GetFrag(n.SeqTwoBit, refEnd-basesToTake, refEnd)...)
}*/

func LeftAlignTraversal(n *Node, seq []dna.Base, refEnd int, currentPath []uint32, ext int, read []dna.Base, scores [][]int64, matrix *MatrixAln, sk scoreKeeper, dynamicScore dynamicScoreKeeper, pool *sync.Pool) ([]cigar.ByteCigar, int64, int, int, []uint32) {
	if len(seq) >= ext {
		log.Fatalf("Error: left traversal, the length=%d of DNA sequence in previous nodes should not be enough to satisfy the desired extenion=%d.\n", len(seq), ext)
	}
	s := pool.Get().(*dnaPool)
	s.Seq = s.Seq[:0]
	s.Seq = getLeftTargetBases(n, ext, refEnd, seq, s.Seq)
	s.Path = s.Path[:0]
	s.Path = append(s.Path, currentPath...)
	AddPath(s.Path, n.Id)
	if len(seq)+refEnd >= ext || len(n.Prev) == 0 {
		sk.leftScore, sk.leftAlignment, sk.minTarget, _, sk.minQuery, _ = LeftDynamicAln(s.Seq, read, scores, matrix, -600, dynamicScore)
		sk.minTarget = refEnd - len(s.Seq) - len(seq) + sk.minTarget
		sk.leftPath = s.Path
		pool.Put(s)
		return sk.leftAlignment, sk.leftScore, sk.minTarget, sk.minQuery, sk.leftPath
	} else {
		//A very negative number
		sk.leftScore = math.MinInt64
		var score int64
		for _, i := range n.Prev {
			dynamicScore.route, score, s.refStart, s.queryStart, s.Path = LeftAlignTraversal(i.Dest, s.Seq, len(i.Dest.Seq), s.Path, ext, read, scores, matrix, sk, dynamicScore, pool)
			if score > sk.leftScore {
				sk.leftScore = score
				sk.leftAlignment = dynamicScore.route
				sk.minTarget = refEnd - len(s.Seq) - len(seq) + s.refStart
				sk.minQuery = s.queryStart
				sk.leftPath = s.Path
			}
		}
		pool.Put(s)
		cigar.ReverseBytesCigar(sk.leftAlignment)
		ReversePath(sk.leftPath)
		return sk.leftAlignment, sk.leftScore, sk.minTarget, sk.minQuery, sk.leftPath
	}
}

func ReversePath(alpha []uint32) {
	var i, off int
	for i, off = len(alpha)/2-1, len(alpha)-1-i; i >= 0; i-- {
		alpha[i], alpha[off] = alpha[off], alpha[i]
	}
}
