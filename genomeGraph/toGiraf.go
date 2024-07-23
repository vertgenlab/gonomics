package genomeGraph

import (
	"fmt"
	"math"
	"strings"
	"sync"

	"github.com/vertgenlab/gonomics/cigar"
	"github.com/vertgenlab/gonomics/fastq"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/giraf"
	"github.com/vertgenlab/gonomics/numbers/parse"
	"github.com/vertgenlab/gonomics/sam"
)

func GraphSmithWatermanToGiraf(gg *GenomeGraph, read fastq.FastqBig, seedHash map[uint64][]uint64, settings *GraphSettings, memoryPool *sync.Pool, sk scoreKeeper, seedBuildHelper *sync.Pool) *giraf.Giraf {
	var currBest giraf.Giraf = giraf.Giraf{
		QName:     read.Name,
		QStart:    0,
		QEnd:      0,
		PosStrand: true,
		Path:      giraf.Path{},
		Cigar:     make([]cigar.Cigar, 0, 1),
		AlnScore:  0,
		MapQ:      255,
		Seq:       read.Seq,
		Qual:      read.Qual,
		Notes:     []giraf.Note{{Tag: []byte{'X', 'O'}, Type: 'Z', Value: "~"}},
	}
	resetScoreKeeper(sk)
	sk.perfectScore = perfectMatchBig(read, settings.ScoreMatrix)
	sk.queryLength = len(read.Seq)
	settings.Extention = int(sk.perfectScore/600) + sk.queryLength

	seedHelper := seedBuildHelper.Get().(*SeedMemory)
	defer seedBuildHelper.Put(seedHelper)

	//settings.MaxMatch, settings.MinMatch, settings.LeastSevereMismatch, settings.LeastSevereMatchMismatchChange = MismatchStats(settings.ScoreMatrix)
	seedHelper.Hits = seedMapMemPool(seedHash, gg.Nodes, &read, sk, settings, seedHelper)

	for i := 0; i < len(seedHelper.Hits); i++ {
		sk.currSeed = &seedHelper.Hits[i]
		sk.tailSeed = getLastPart(sk.currSeed)
		if sk.currSeed.PosStrand {
			sk.currSeq = read.Seq
		} else {
			sk.currSeq = read.SeqRc
		}
		sk.seedScore = scoreSeedSeq(sk.currSeq, sk.currSeed.QueryStart, sk.tailSeed.QueryStart+sk.tailSeed.Length, settings.ScoreMatrix)
		if int(sk.currSeed.TotalLength) == sk.queryLength {
			sk.targetStart = int(sk.currSeed.TargetStart)
			sk.targetEnd = int(sk.tailSeed.TargetStart + sk.tailSeed.Length)
			sk.queryStart = int(sk.currSeed.QueryStart)
			sk.currScore = sk.seedScore
		} else {
			sk.leftAlignment, sk.leftScore, sk.targetStart, sk.queryStart, sk.leftPath = LeftAlignTraversal(&gg.Nodes[sk.currSeed.TargetId], sk.leftSeq, int(sk.currSeed.TargetStart), sk.leftPath, sk.currSeq[:sk.currSeed.QueryStart], settings, sk, memoryPool)
			sk.rightAlignment, sk.rightScore, sk.targetEnd, sk.queryEnd, sk.rightPath = RightAlignTraversal(&gg.Nodes[sk.tailSeed.TargetId], sk.rightSeq, int(sk.tailSeed.TargetStart+sk.tailSeed.Length), sk.rightPath, sk.currSeq[sk.tailSeed.QueryStart+sk.tailSeed.Length:], settings, sk, memoryPool)
			sk.currScore = sk.leftScore + sk.seedScore + sk.rightScore
		}
		if sk.currScore > int64(currBest.AlnScore) {
			currBest.QStart = sk.queryStart
			currBest.QEnd = int(sk.currSeed.QueryStart) + currBest.QStart + sk.queryEnd + int(sk.currSeed.TotalLength) - 1
			currBest.PosStrand = sk.currSeed.PosStrand
			currBest.Path = setPath(currBest.Path, sk.targetStart, CatPaths(CatPaths(sk.leftPath, getSeedPath(sk.currSeed)), sk.rightPath), sk.targetEnd)
			currBest.Cigar = cigar.CatCigar(cigar.AddCigar(sk.leftAlignment, cigar.Cigar{RunLength: int(sk.currSeed.TotalLength), Op: cigar.Match}), sk.rightAlignment)
			currBest.AlnScore = int(sk.currScore)
			currBest.Seq = sk.currSeq
		}
	}
	if len(currBest.Cigar) > 0 {
		currBest.Cigar = cigar.SoftClipBases(currBest.QStart, sk.queryLength, currBest.Cigar)
	}

	if !currBest.PosStrand {
		fastq.ReverseQualUint8Record(currBest.Qual)
	}
	return &currBest
}

func readFastqGsw(fileOne string, fileTwo string, answer chan<- fastq.PairedEndBig) {
	readOne, readTwo := fileio.NewByteReader(fileOne), fileio.NewByteReader(fileTwo)
	for fq, done := fastq.ReadFqBigPair(readOne, readTwo); !done; fq, done = fastq.ReadFqBigPair(readOne, readTwo) {
		answer <- fq
	}
	close(answer)
}

func WrapPairGiraf(gg *GenomeGraph, fq fastq.PairedEndBig, seedHash map[uint64][]uint64, settings *GraphSettings, memoryPool *sync.Pool, sk scoreKeeper, seedBuildHelper *sync.Pool) giraf.GirafPair {
	var mappedPair giraf.GirafPair = giraf.GirafPair{
		Fwd: *GraphSmithWatermanToGiraf(gg, fq.Fwd, seedHash, settings, memoryPool, sk, seedBuildHelper),
		Rev: *GraphSmithWatermanToGiraf(gg, fq.Rev, seedHash, settings, memoryPool, sk, seedBuildHelper),
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
