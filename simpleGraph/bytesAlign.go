package simpleGraph

import(
	"fmt"
	"github.com/vertgenlab/gonomics/common"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/fastq"
	"github.com/vertgenlab/gonomics/giraf"
	"log"
	"strings"
	"sync"
	"github.com/edotau/simpleio"
)

func SimplyGsw(gg *SimpleGraph, read *fastq.FastqBig, seedHash map[uint64][]uint64, seedLen int, stepSize int, scoreMatrix [][]int64, m [][]int64, trace [][]simpleio.CigarOp) *giraf.Giraf {
	var currBest giraf.Giraf = giraf.Giraf{
		QName:     read.Name,
		QStart:    0,
		QEnd:      0,
		PosStrand: true,
		Path:      &giraf.Path{},
		Aln:       []simpleio.ByteCigar{},
		AlnScore:  0,
		MapQ:      255,
		Seq:       read.Seq,
		Qual:      read.Qual,
		Notes:     []giraf.Note{giraf.Note{Tag: "XO", Type: 'Z', Value: "~"}},
	}
	var leftAlignment, rightAlignment []simpleio.ByteCigar = []simpleio.ByteCigar{}, []simpleio.ByteCigar{}
	var minTarget, maxTarget int
	var minQuery, maxQuery int
	var leftScore, rightScore int64 = 0, 0
	var leftPath, rightPath []uint32
	var currScore int64 = 0
	perfectScore := perfectMatchBig(read, scoreMatrix)
	extension := int(perfectScore/600) + len(read.Seq)
	var seeds []*SeedDev
	seeds = findSeedsInSmallMapWithMemPool(seedHash, gg.Nodes, read, seedLen, perfectScore, scoreMatrix)
	SortSeedDevByLen(seeds)
	var tailSeed *SeedDev
	var seedScore int64
	var currSeq []dna.Base
	var currSeed *SeedDev

	for i := 0; i < len(seeds) && seedCouldBeBetter(int64(seeds[i].TotalLength), int64(currBest.AlnScore), perfectScore, int64(len(read.Seq)), 100, 90, -196, -296); i++ {
		currSeed = seeds[i]
		tailSeed = getLastPart(currSeed)
		if currSeed.PosStrand {
			currSeq = read.Seq
		} else {
			currSeq = read.SeqRc
		}
		seedScore = scoreSeedSeq(currSeq, currSeed.QueryStart, tailSeed.QueryStart+tailSeed.Length, scoreMatrix)
		if int(currSeed.TotalLength) == len(currSeq) {
			currScore = seedScore
			minTarget = int(currSeed.TargetStart)
			maxTarget = int(tailSeed.TargetStart + tailSeed.Length)
			minQuery = int(currSeed.QueryStart)
			maxQuery = int(currSeed.TotalLength - 1)
		} else {
			leftAlignment, leftScore, minTarget, minQuery, leftPath = LeftAlignTraversal(gg.Nodes[currSeed.TargetId], []dna.Base{}, int(currSeed.TargetStart), []uint32{}, extension-int(currSeed.TotalLength), currSeq[:currSeed.QueryStart], m, trace)
			rightAlignment, rightScore, maxTarget, maxQuery, rightPath = RightAlignTraversal(gg.Nodes[tailSeed.TargetId], []dna.Base{}, int(tailSeed.TargetStart+tailSeed.Length), []uint32{}, extension-int(currSeed.TotalLength), currSeq[tailSeed.QueryStart+tailSeed.Length:], m, trace)
			//log.Printf("left alignment: %s\n", simpleio.ByteCigarString(leftAlignment))
			simpleio.ByteCigarString(rightAlignment)
			//log.Printf("right alignment: %s\n", simpleio.ByteCigarString(rightAlignment))
		}
		currScore = leftScore + seedScore + rightScore
		if currScore > int64(currBest.AlnScore) {
			currBest.QStart = minQuery
			currBest.QEnd = maxQuery
			currBest.PosStrand = currSeed.PosStrand
			currBest.Path = setPath(currBest.Path, minTarget, CatPaths(CatPaths(leftPath, getSeedPath(currSeed)), rightPath), maxTarget)
			//currBest.Aln = AddSClip(minQuery, len(currSeq), cigar.CatCigar(cigar.AddCigar(leftAlignment, &cigar.Cigar{RunLength: int64(sumLen(currSeed)), Op: 'M'}), rightAlignment))
			
			currBest.Aln = simpleio.SoftClipBases(minQuery, len(currSeq), simpleio.CatByteCigar(simpleio.AddCigar(leftAlignment, simpleio.ByteCigar{RunLen: sumLen(currSeed), Op: 'M'}), rightAlignment))
			currBest.AlnScore = int(currScore)
			currBest.Seq = currSeq
			if gg.Nodes[currBest.Path.Nodes[0]].Info != nil {
				currBest.Notes[0].Value = fmt.Sprintf("%s=%d", gg.Nodes[currBest.Path.Nodes[0]].Name, gg.Nodes[currBest.Path.Nodes[0]].Info.Start)
				currBest.Notes = append(currBest.Notes, infoToNotes(gg.Nodes, currBest.Path.Nodes))
			} else {
				currBest.Notes[0].Value = fmt.Sprintf("%s=%d", gg.Nodes[currBest.Path.Nodes[0]].Name, 1)
			}
		}
	}
	if !currBest.PosStrand {
		fastq.ReverseQualUint8Record(currBest.Qual)
	}
	return &currBest
}


func RightAlignTraversal(n *Node, seq []dna.Base, start int, currentPath []uint32, ext int, read []dna.Base, m [][]int64, trace [][]simpleio.CigarOp) ([]simpleio.ByteCigar, int64, int, int, []uint32) {
	currentPath = append(currentPath, n.Id)
	var bestTargetEnd, bestQueryEnd, targetEnd, queryEnd int
	var bestScore, score int64
	var bestAlignment, alignment []simpleio.ByteCigar
	var path, bestPath []uint32

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
		return alignment, score, targetEnd, queryEnd, currentPath
	} else {
		bestScore = -1
		tmpPath := make([]uint32, len(currentPath))
		copy(tmpPath, currentPath)
		for _, i := range n.Next {
			AddPath(i.Dest.Id, currentPath)
			alignment, score, targetEnd, queryEnd, path = RightAlignTraversal(i.Dest, s, 0, currentPath, ext, read, m, trace)
			if score > bestScore {
				bestScore = score
				bestAlignment = alignment
				bestTargetEnd = targetEnd
				bestQueryEnd = queryEnd
				bestPath = path
			}
		}
	}
	return bestAlignment, bestScore, bestTargetEnd, bestQueryEnd, bestPath
}

func LeftAlignTraversal(n *Node, seq []dna.Base, refEnd int, currentPath []uint32, ext int, read []dna.Base, m [][]int64, trace [][]simpleio.CigarOp) ([]simpleio.ByteCigar, int64, int, int, []uint32) {
	currentPath = append(currentPath, n.Id)
	var bestQueryStart, queryStart, refStart, bestRefStart int
	var bestScore, score int64
	var bestAlignment, alignment []simpleio.ByteCigar
	var path, bestPath []uint32

	var availableBases int = len(seq) + refEnd
	var targetLength int = common.Min(availableBases, ext)
	var basesToTake int = targetLength - len(seq)
	var s []dna.Base = make([]dna.Base, targetLength)
	copy(s[0:basesToTake], n.Seq[refEnd-basesToTake:refEnd])
	copy(s[basesToTake:targetLength], seq)
	if availableBases >= ext || len(n.Next) == 0 {
		score, alignment, refStart, _, queryStart, _ = simpleio.LeftDynamicAln(s, read, HumanChimpTwoScoreMatrix, -600, m, trace)
		reversePath(currentPath)
		return alignment, score, refEnd - basesToTake + refStart, queryStart, currentPath
	} else {
		bestScore = -1
		for _, i := range n.Prev {
			tmp := make([]uint32, len(currentPath))
			copy(tmp, currentPath)
			alignment, score, refStart, queryStart, path = LeftAlignTraversal(i.Dest, s, len(i.Dest.Seq), currentPath, ext, read, m, trace)
			if score > bestScore {
				bestScore = score
				bestAlignment = alignment
				bestRefStart = refStart
				bestQueryStart = queryStart
				bestPath = path
			}
		}
	}
	return bestAlignment, bestScore, bestRefStart, bestQueryStart, bestPath
}

func RightDynamicAln(alpha []dna.Base, beta []dna.Base, scores [][]int64, gapPen int64, m [][]int64, trace [][]simpleio.CigarOp) (int64, []simpleio.ByteCigar, int, int, int, int) {
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

func LeftDynamicAln(alpha []dna.Base, beta []dna.Base, scores [][]int64, gapPen int64, m [][]int64, trace [][]simpleio.CigarOp) (int64, []simpleio.ByteCigar, int, int, int, int) {
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

func RoutineSimplyGiraf(gg *SimpleGraph, seedHash map[uint64][]uint64, seedLen int, stepSize int, scoreMatrix [][]int64, inputChan <-chan *fastq.PairedEndBig, outputChan chan<- *giraf.GirafPair, wg *sync.WaitGroup) {
	m, trace := simpleio.MatrixSetup(10000)
	for read := range inputChan {
		outputChan <- WrapSimplyGsw(gg, read, seedHash, seedLen, stepSize, scoreMatrix, m, trace)
	}
	wg.Done()
}

func WrapSimplyGsw(gg *SimpleGraph, readPair *fastq.PairedEndBig, seedHash map[uint64][]uint64, seedLen int, stepSize int, scoreMatrix [][]int64, m [][]int64, trace [][]simpleio.CigarOp) *giraf.GirafPair {
	var mappedPair giraf.GirafPair = giraf.GirafPair{Fwd: nil, Rev: nil}
	mappedPair.Fwd = SimplyGsw(gg, readPair.Fwd, seedHash, seedLen, stepSize, scoreMatrix, m, trace)
	mappedPair.Rev = SimplyGsw(gg, readPair.Rev, seedHash, seedLen, stepSize, scoreMatrix, m, trace)
	setGirafFlags(&mappedPair)
	return &mappedPair
}

func SoftClipBases(front int, lengthOfRead int, cig []simpleio.ByteCigar) []simpleio.ByteCigar {
	var runLen int = simpleio.QueryRunLen(cig)
	if runLen < lengthOfRead {
		answer := make([]simpleio.ByteCigar, 0, len(cig)+2)
		if front > 0 {
			answer = append(answer, simpleio.ByteCigar{RunLen: uint32(front), Op: 'S'})
		}
		answer = append(answer, cig...)
		if front+simpleio.QueryRunLen(cig) < lengthOfRead {
			answer = append(answer, simpleio.ByteCigar{RunLen: uint32(lengthOfRead-front - runLen), Op: 'S'})
		}
		return answer
	} else {
		return cig
	}
}

func IsGirafCorrect(aln *giraf.Giraf, genome *SimpleGraph) bool {
	qName := strings.Split(aln.QName, "_")
	if len(qName) < 6 {
		log.Fatalf("Error: input giraf file does not match simulation format...\n")
	}
	if common.StringToInt(qName[1]) == aln.Path.TStart && common.StringToUint32(qName[1]) == aln.Path.Nodes[0] {
		return true
	} else {
		return false
	}
}
