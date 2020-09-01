package simpleGraph

import (
	"fmt"
	"github.com/vertgenlab/gonomics/cigar"
	"github.com/vertgenlab/gonomics/common"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/fastq"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/giraf"
	"github.com/vertgenlab/gonomics/sam"
	"math"
	"strings"
	"sync"
)

func GraphSmithWatermanToGiraf(gg *SimpleGraph, read fastq.FastqBig, seedHash map[uint64][]uint64, seedLen int, stepSize int, matrix *MatrixAln, scoreMatrix [][]int64, seedPool *sync.Pool, dnaPool *sync.Pool, sk scoreKeeper, dynamicScore dynamicScoreKeeper) *giraf.Giraf {
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
	return &currBest
}

func readFastqGsw(fileOne string, fileTwo string, answer chan<- fastq.PairedEndBig) {
	readOne, readTwo := fileio.NewSimpleReader(fileOne), fileio.NewSimpleReader(fileTwo)
	for fq, done := fastq.ReadFqBigPair(readOne, readTwo); !done; fq, done = fastq.ReadFqBigPair(readOne, readTwo) {
		answer <- *fq
	}
	close(answer)
}

/*
func GirafToExplicitCigar(giraf *giraf.Giraf, graph *SimpleGraph) []*cigar.Cigar {
	var answer []*cigar.Cigar
	var seqIdx, refIdx, pathIdx int
	refIdx = giraf.Path.TStart
	var k, runLenCount int64

	if giraf.Aln[0].Op == '*' {
		return nil
	}

	for i := 0; i < len(giraf.Aln); i++ {
		switch giraf.Aln[i].Op {
		case 'M':
			runLenCount = 0

			for k = 0; k < giraf.Aln[i].RunLength; k++ {
				if refIdx > len(graph.Nodes[giraf.Path.Nodes[pathIdx]].Seq)-1 {
					pathIdx++
					refIdx = 0
				}
				if giraf.Seq[seqIdx] == graph.Nodes[giraf.Path.Nodes[pathIdx]].Seq[refIdx] {
					runLenCount++
				} else {
					if runLenCount > 0 {
						// Append the matching bases so far
						answer = append(answer, &cigar.Cigar{RunLength: runLenCount, Op: '='})
					}
					// Append the mismatch base
					if answer == nil {
						answer = append(answer, &cigar.Cigar{RunLength: 1, Op: 'X', Sequence: []dna.Base{giraf.Seq[k]}})
					} else if answer[len(answer)-1].Op == 'X' {
						answer[len(answer)-1].RunLength++
						answer[len(answer)-1].Sequence = append(answer[len(answer)-1].Sequence, giraf.Seq[k])
					} else {
						answer = append(answer, &cigar.Cigar{RunLength: 1, Op: 'X', Sequence: []dna.Base{giraf.Seq[k]}})
					}
					runLenCount = 0
				}
				seqIdx++
				refIdx++
			}

			if runLenCount > 0 {
				answer = append(answer, &cigar.Cigar{RunLength: runLenCount, Op: '='})
			}

		case 'I':
			var insSeq []dna.Base
			for k = 0; k < giraf.Aln[i].RunLength; k++ {
				insSeq = append(insSeq, giraf.Seq[seqIdx])
				seqIdx++
			}
			answer = append(answer, &cigar.Cigar{RunLength: giraf.Aln[i].RunLength, Op: 'I', Sequence: insSeq})

		case 'X':
			log.Println("WARNING: The input cigar already has explicit formatting")
			return giraf.Aln

		case '=':
			log.Println("WARNING: The input cigar already has explicit formatting")
			return giraf.Aln

		default:
			answer = append(answer, giraf.Aln[i])
			if cigar.ConsumesReference(giraf.Aln[i].Op) {
				refIdx += int(giraf.Aln[i].RunLength)
			}
			if cigar.ConsumesQuery(giraf.Aln[i].Op) {
				seqIdx += int(giraf.Aln[i].RunLength)
			}
		}
	}
	return answer
}*/

type ScoreMatrixHelper struct {
	Matrix                         [][]int64
	MaxMatch                       int64
	MinMatch                       int64
	LeastSevereMismatch            int64
	LeastSevereMatchMismatchChange int64
}

func getScoreMatrixHelp(scoreMatrix [][]int64) *ScoreMatrixHelper {
	help := ScoreMatrixHelper{Matrix: scoreMatrix}
	help.MaxMatch, help.MinMatch, help.LeastSevereMismatch, help.LeastSevereMatchMismatchChange = MismatchStats(scoreMatrix)
	return &help
}

func MismatchStats(scoreMatrix [][]int64) (int64, int64, int64, int64) {
	var maxMatch int64 = 0
	var minMatch int64
	var leastSevereMismatch int64 = scoreMatrix[0][1]
	var i, j int
	for i = 0; i < len(scoreMatrix); i++ {
		for j = 0; j < len(scoreMatrix[i]); j++ {
			if scoreMatrix[i][j] > maxMatch {
				minMatch = maxMatch
				maxMatch = scoreMatrix[i][j]
			} else {
				if scoreMatrix[i][j] < 0 && leastSevereMismatch < scoreMatrix[i][j] {
					leastSevereMismatch = scoreMatrix[i][j]
				}
			}

		}
	}
	var leastSevereMatchMismatchChange int64 = leastSevereMismatch - maxMatch
	return maxMatch, minMatch, leastSevereMismatch, leastSevereMatchMismatchChange
}

func WrapPairGiraf(gg *SimpleGraph, fq fastq.PairedEndBig, seedHash map[uint64][]uint64, seedLen int, stepSize int, matrix *MatrixAln, scoreMatrix [][]int64, seedPool *sync.Pool, dnaPool *sync.Pool, sk scoreKeeper, dynamicScore dynamicScoreKeeper) giraf.GirafPair {
	var mappedPair giraf.GirafPair = giraf.GirafPair{
		Fwd: *GraphSmithWatermanToGiraf(gg, fq.Fwd, seedHash, seedLen, stepSize, matrix, scoreMatrix, seedPool, dnaPool, sk, dynamicScore),
		Rev: *GraphSmithWatermanToGiraf(gg, fq.Rev, seedHash, seedLen, stepSize, matrix, scoreMatrix, seedPool, dnaPool, sk, dynamicScore),
	}
	//setGirafFlags(&mappedPair)
	return mappedPair
}

// setGirafFlags generates the appropriate flags for each giraf in a pair
func setGirafFlags(pair *giraf.GirafPair) {
	pair.Fwd.Flag = getGirafFlags(&pair.Fwd)
	pair.Rev.Flag = getGirafFlags(&pair.Rev)
	pair.Fwd.Flag += 8  // Forward
	pair.Fwd.Flag += 16 // Paired Reads
	pair.Fwd.Flag += 16 // Paired Reads
	if isProperPairAlign(pair) {
		pair.Fwd.Flag += 1 // Properly Aligned
		pair.Rev.Flag += 1 // Properly Aligned
	}
}

func GirafToSam(ag *giraf.Giraf) *sam.SamAln {
	curr := sam.SamAln{QName: ag.QName, Flag: 4, RName: "*", Pos: 0, MapQ: 255, Cigar: []*cigar.Cigar{&cigar.Cigar{Op: '*'}}, RNext: "*", PNext: 0, TLen: 0, Seq: ag.Seq, Qual: fastq.QualString(ag.Qual), Extra: "BZ:i:0\tGP:Z:-1\tXO:Z:~"}
	//read is unMapped
	if strings.Compare(ag.Notes[0].Value, "~") == 0 {
		return &curr
	} else {
		target := strings.Split(ag.Notes[0].Value, "=")
		curr.RName = target[0]
		curr.Pos = int64(ag.Path.TStart) + common.StringToInt64(target[1])
		curr.Flag = getSamFlags(ag)
		if len(ag.Notes) == 2 {
			curr.Extra = fmt.Sprintf("BZ:i:%d\tGP:Z:%s\tXO:i:%d\t%s", ag.AlnScore, PathToString(ag.Path.Nodes), ag.Path.TStart, giraf.NoteToString(ag.Notes[1]))
		} else {
			curr.Extra = fmt.Sprintf("BZ:i:%d\tGP:Z:%s\tXO:i:%d", ag.AlnScore, PathToString(ag.Path.Nodes), ag.Path.TStart)
		}
	}
	return &curr
}

func GirafPairToSam(ag giraf.GirafPair) *sam.PairedSamAln {
	var mappedPair sam.PairedSamAln = sam.PairedSamAln{FwdSam: &sam.SamAln{}, RevSam: &sam.SamAln{}}
	mappedPair.FwdSam = GirafToSam(&ag.Fwd)
	mappedPair.RevSam = GirafToSam(&ag.Rev)
	mappedPair.FwdSam.Flag += 64 + 2
	mappedPair.RevSam.Flag += 128 + 2
	if isProperPairAlign(&ag) {
		mappedPair.FwdSam.Flag += 1
		mappedPair.RevSam.Flag += 1
	}
	return &mappedPair
}

func isProperPairAlign(mappedPair *giraf.GirafPair) bool {
	if math.Abs(float64(mappedPair.Fwd.Path.TStart-mappedPair.Rev.Path.TStart)) < 10000 {
		if mappedPair.Fwd.Path.TStart < mappedPair.Rev.Path.TStart && mappedPair.Fwd.PosStrand && !mappedPair.Rev.PosStrand {
			return true
		}
		if mappedPair.Fwd.Path.TStart > mappedPair.Rev.Path.TStart && !mappedPair.Fwd.PosStrand && mappedPair.Rev.PosStrand {
			return true
		}
	}
	return false
}

func getGirafFlags(ag *giraf.Giraf) uint8 {
	var answer uint8
	if ag.PosStrand {
		answer += 4 // Positive Strand
	}
	if ag.AlnScore < 1200 {
		answer += 2 // Unmapped
	}
	return answer
}

func getSamFlags(ag *giraf.Giraf) int64 {
	var answer int64
	if !ag.PosStrand {
		answer += 16
	}
	if ag.AlnScore < 1200 {
		answer += 4
	}
	return answer
}

func setPath(p giraf.Path, targetStart int, nodes []uint32, targetEnd int) giraf.Path {
	p.TStart = targetStart
	p.Nodes = nodes
	p.TEnd = targetEnd
	return p
}

func vInfoToValue(n *Node) string {
	var answer string
	switch {
	case n.Info.Variant == 1:
		answer = fmt.Sprintf("%d=%s", n.Id, "snp")
	case n.Info.Variant == 2:
		answer = fmt.Sprintf("%d=%s", n.Id, "ins")
	case n.Info.Variant == 3:
		answer = fmt.Sprintf("%d=%s", n.Id, "del")
	}
	return answer
}

func infoToNotes(nodes []*Node, path []uint32) giraf.Note {
	var vInfo giraf.Note = giraf.Note{Tag: []byte{'X', 'V'}, Type: 'Z'}
	vInfo.Value = fmt.Sprintf("%d_%d", nodes[0].Info.Allele, nodes[0].Info.Variant)
	if len(path) > 0 {
		for i := 1; i < len(path); i++ {
			if nodes[i].Info.Variant > 0 {
				vInfo.Value += fmt.Sprintf(",%s", vInfoToValue(nodes[path[i]]))
			} else {
				vInfo.Value += fmt.Sprintf(",%d_%d", nodes[i].Info.Allele, nodes[path[i]].Info.Variant)
			}

		}
	}
	return vInfo
}
