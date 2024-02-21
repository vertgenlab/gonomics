package genomeGraph

import (
	"fmt"
	"math"
	"strings"
	"sync"

	"github.com/vertgenlab/gonomics/cigar"
	"github.com/vertgenlab/gonomics/fastq"
	"github.com/vertgenlab/gonomics/giraf"
	"github.com/vertgenlab/gonomics/numbers/parse"
	"github.com/vertgenlab/gonomics/sam"
)

func GraphSmithWatermanToGiraf(gg *GenomeGraph, read fastq.FastqBig, seedHash map[uint64][]uint64, seedLen int, stepSize int, matrix *MatrixAln, scoreMatrix [][]int64, seedPool *sync.Pool, dnaPool *sync.Pool, sk scoreKeeper, dynamicScore dynamicScoreKeeper, seedBuildHelper *seedHelper) *giraf.Giraf {
	seeds := seedPool.Get().(*memoryPool)
	defer seedPool.Put(seeds)

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
		Notes:     []giraf.Note{{Tag: []byte{cigar.Mismatch, 'O'}, Type: 'Z', Value: "~"}},
	}

	sk.perfectScore = perfectMatchBig(read, scoreMatrix)
	sk.extension = int(sk.perfectScore/600) + len(read.Seq)

	seeds.Hits = seeds.Hits[:0]
	seeds.Worker = seeds.Worker[:0]
	seeds.Hits = seedMapMemPool(seedHash, gg.Nodes, &read, seedLen, sk.perfectScore, scoreMatrix, seeds.Hits, seeds.Worker, seedBuildHelper)

	for i := 0; i < len(seeds.Hits) && seedCouldBeBetter(int64(seeds.Hits[i].TotalLength), int64(currBest.AlnScore), sk.perfectScore, int64(len(read.Seq)), 100, 90, -196, -296); i++ {
		sk.currSeed = seeds.Hits[i]
		sk.tailSeed = *getLastPart(&sk.currSeed)
		if sk.currSeed.PosStrand {
			sk.currSeq = read.Seq
		} else {
			sk.currSeq = read.SeqRc
		}
		sk.seedScore = scoreSeedSeq(sk.currSeq, sk.currSeed.QueryStart, sk.tailSeed.QueryStart+sk.tailSeed.Length, scoreMatrix)
		if int(sk.currSeed.TotalLength) == len(sk.currSeq) {
			sk.targetStart = int(sk.currSeed.TargetStart)
			sk.targetEnd = int(sk.tailSeed.TargetStart + sk.tailSeed.Length)
			sk.queryStart = int(sk.currSeed.QueryStart)
			sk.currScore = sk.seedScore
		} else {
			sk.leftAlignment, sk.leftScore, sk.targetStart, sk.queryStart, sk.leftPath = AlignGraphTraversal(&gg.Nodes[sk.currSeed.TargetId], sk.leftSeq, int(sk.currSeed.TargetStart), sk.leftPath, sk.extension-int(sk.currSeed.TotalLength), sk.currSeq[:sk.currSeed.QueryStart], scoreMatrix, matrix, &sk, &dynamicScore, dnaPool, leftTraversal)
			sk.rightAlignment, sk.rightScore, sk.targetEnd, sk.queryEnd, sk.rightPath = AlignGraphTraversal(&gg.Nodes[sk.tailSeed.TargetId], sk.rightSeq, int(sk.tailSeed.TargetStart+sk.tailSeed.Length), sk.rightPath, sk.extension-int(sk.currSeed.TotalLength), sk.currSeq[sk.tailSeed.QueryStart+sk.tailSeed.Length:], scoreMatrix, matrix, &sk, &dynamicScore, dnaPool, rightTraversal)
			sk.currScore = sk.leftScore + sk.seedScore + sk.rightScore
		}
		if sk.currScore > int64(currBest.AlnScore) {
			currBest.QStart = sk.queryStart
			currBest.QEnd = int(sk.currSeed.QueryStart) + sk.queryStart + sk.queryEnd + int(sk.currSeed.TotalLength) - 1
			currBest.PosStrand = sk.currSeed.PosStrand
			currBest.Path = setPath(currBest.Path, sk.targetStart, CatPaths(CatPaths(sk.leftPath, getSeedPath(&sk.currSeed)), sk.rightPath), sk.targetEnd)
			currBest.Cigar = cigar.AppendSoftClipBases(sk.queryStart, len(sk.currSeq), cigar.CatByteCigar(cigar.AddCigarByte(sk.leftAlignment, cigar.ByteCigar{RunLen: uint16(sk.currSeed.TotalLength), Op: 'M'}), sk.rightAlignment))
			currBest.AlnScore = int(sk.currScore)
			currBest.Seq = sk.currSeq
		}
	}

	if !currBest.PosStrand {
		fastq.ReverseQualUint8Record(currBest.Qual)
	}
	return &currBest
}

type ScoreMatrixHelper struct {
	Matrix                         [][]int64
	MaxMatch                       int64
	MinMatch                       int64
	LeastSevereMismatch            int64
	LeastSevereMatchMismatchChange int64
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

func WrapPairGiraf(gg *GenomeGraph, fq fastq.PairedEndBig, seedHash map[uint64][]uint64, seedLen int, stepSize int, matrix *MatrixAln, scoreMatrix [][]int64, seedPool *sync.Pool, dnaPool *sync.Pool, sk scoreKeeper, dynamicScore dynamicScoreKeeper, seedBuildHelper *seedHelper) giraf.GirafPair {
	var mappedPair giraf.GirafPair = giraf.GirafPair{
		Fwd: *GraphSmithWatermanToGiraf(gg, fq.Fwd, seedHash, seedLen, stepSize, matrix, scoreMatrix, seedPool, dnaPool, sk, dynamicScore, seedBuildHelper),
		Rev: *GraphSmithWatermanToGiraf(gg, fq.Rev, seedHash, seedLen, stepSize, matrix, scoreMatrix, seedPool, dnaPool, sk, dynamicScore, seedBuildHelper),
	}
	setGirafFlags(&mappedPair)
	return mappedPair
}

// setGirafFlags generates the appropriate flags for each giraf in a pair.
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

func GirafToSam(ag *giraf.Giraf) sam.Sam {
	curr := sam.Sam{QName: ag.QName, Flag: 4, RName: "*", Pos: 0, MapQ: 255, Cigar: []cigar.Cigar{{Op: '*'}}, RNext: "*", PNext: 0, TLen: 0, Seq: ag.Seq, Qual: fastq.QualString(ag.Qual), Extra: "BZ:i:0\tGP:Z:-1\tXO:Z:~"}
	//read is unMapped
	if strings.Compare(ag.Notes[0].Value, "~") == 0 {
		return curr
	} else {
		target := strings.Split(ag.Notes[0].Value, "=")
		curr.RName = target[0]
		curr.Pos = uint32(ag.Path.TStart + parse.StringToInt(target[1]))
		curr.Flag = getSamFlags(ag)
		if len(ag.Notes) == 2 {
			curr.Extra = fmt.Sprintf("BZ:i:%d\tGP:Z:%s\tXO:i:%d\t%s", ag.AlnScore, PathToString(ag.Path.Nodes), ag.Path.TStart, giraf.NoteToString(ag.Notes[1]))
		} else {
			curr.Extra = fmt.Sprintf("BZ:i:%d\tGP:Z:%s\tXO:i:%d", ag.AlnScore, PathToString(ag.Path.Nodes), ag.Path.TStart)
		}
	}
	return curr
}

func GirafPairToSam(ag giraf.GirafPair) sam.MatePair {
	var mappedPair sam.MatePair = sam.MatePair{Fwd: sam.Sam{}, Rev: sam.Sam{}}
	mappedPair.Fwd = GirafToSam(&ag.Fwd)
	mappedPair.Rev = GirafToSam(&ag.Rev)
	mappedPair.Fwd.Flag += 64 + 2
	mappedPair.Rev.Flag += 128 + 2
	if isProperPairAlign(&ag) {
		mappedPair.Fwd.Flag += 1
		mappedPair.Rev.Flag += 1
	}
	return mappedPair
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

func getSamFlags(ag *giraf.Giraf) uint16 {
	var answer uint16
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

/*func vInfoToValue(n *Node) string {
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
}*/

/*func infoToNotes(nodes []*Node, path []uint32) giraf.Note {
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
}*/
