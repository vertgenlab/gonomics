package chain

import (
	"github.com/vertgenlab/gonomics/axt"
	"github.com/vertgenlab/gonomics/common"
	"github.com/vertgenlab/gonomics/dna"
	"log"
)
//TODO: Next to figure out a better way to pass in target and query
func ChainToAxt(ch *Chain, target []dna.Base, query []dna.Base) *axt.Axt {
	var tLen, qLen int = ch.TEnd - ch.TStart, ch.QEnd - ch.QStart
	var answer *axt.Axt = &axt.Axt{
		RName:      ch.TName,
		RStart:     int64(ch.TStart) + 1,
		REnd:       int64(ch.TEnd) + 1,
		QName:      ch.QName,
		QStart:     int64(ch.QStart) + 1,
		QEnd:       int64(ch.QEnd) + 1,
		QStrandPos: ch.QStrand,
		Score:      int64(ch.Score),
		RSeq:       make([]dna.Base, 0, tLen),
		QSeq:       make([]dna.Base, 0, qLen),
	}
	var tIndex, qIndex int = ch.TStart, ch.QStart
	for _, each := range ch.Alignment {
		answer.RSeq = append(answer.RSeq, getSequence(target, tIndex, each.Size)...)
		answer.QSeq = append(answer.QSeq, getSequence(query, qIndex, each.Size)...)

		tIndex, qIndex = tIndex+each.Size, qIndex+each.Size
		//update idx
		if each.TBases > 0 {
			answer.RSeq = append(answer.RSeq, getSequence(target, tIndex, each.TBases)...)
			answer.QSeq = append(answer.QSeq, dna.CreateAllGaps(int64(each.TBases))...)
			tIndex += each.TBases
		}
		if each.QBases > 0 {
			answer.QSeq = append(answer.QSeq, getSequence(query, qIndex, each.QBases)...)
			answer.RSeq = append(answer.RSeq, dna.CreateAllGaps(int64(each.QBases))...)
			qIndex += each.QBases
		}
	}
	if len(answer.RSeq) != tLen || len(answer.QSeq) != qLen {
		log.Fatalf("Error: %s indices for target and/or query, may be off by 1\n", axt.ToString(answer, 0))
	}
	return answer
}

//helper function to quickly get sequence at a starting pos plus length this way i dont have to keep doing start:start+length everytime
func getSequence(seq []dna.Base, start int, length int) []dna.Base {
	return seq[start : start+length]
}

//get alignment blocks for the 3 columns containing alignment data
//length of of target and query seqs in axts should be the same i believe, so we can loop over either one.
func getChainCounts(rSeq []dna.Base, qSeq []dna.Base) int {
	var answer int
	for i := 0; i < len(rSeq); i++ {
		if rSeq[i] != dna.Gap && qSeq[i] != dna.Gap {
			answer++
		} else {
			break
		}
	}
	return answer
}

func calcMissingBases(rSeq []dna.Base, qSeq []dna.Base) (int, int) {
	var target, query int
	for i, j := 0, 0; i < len(rSeq) && j < len(qSeq) && (rSeq[i] == dna.Gap || qSeq[i] == dna.Gap); {
		if rSeq[i] == dna.Gap {
			query++
		}
		if qSeq[i] == dna.Gap {
			target++
		}
		i++
		j++
	}
	return target, query
}

func CalcEntireBlock(rSeq []dna.Base, qSeq []dna.Base) []*BaseStats {
	var answer []*BaseStats
	var curr *BaseStats
	if len(rSeq) != len(qSeq) {
		log.Fatalf("Error input sequences should match in length\n")
	}
	for i := 0; i < len(rSeq) && i < len(qSeq); {
		//first count matching seq
		curr = &BaseStats{
			Size:   getChainCounts(rSeq[i:], qSeq[i:]),
			TBases: 0,
			QBases: 0,
		}
		i += curr.Size
		if i == len(rSeq) {
			answer = append(answer, curr)
		} else {
			curr.TBases, curr.QBases = calcMissingBases(rSeq[i:], qSeq[i:])
			answer = append(answer, curr)
			i += common.Max(curr.TBases, curr.QBases)
		}
	}
	return answer
}

func AxtToChain(align *axt.Axt, id int) *Chain {
	var tLen, qLen int = int(align.REnd - align.RStart), int(align.QEnd - align.QStart)
	var answer *Chain = &Chain{
		Score:     int(align.Score),
		TName:     align.RName,
		TSize:     tLen,
		TStrand:   true,
		TStart:    int(align.RStart) - 1,
		TEnd:      int(align.REnd) - 1,
		QName:     align.QName,
		QSize:     qLen,
		QStrand:   align.QStrandPos,
		QStart:    int(align.QStart) - 1,
		QEnd:      int(align.QEnd) - 1,
		Alignment: make([]*BaseStats, 0),
		Id:        id,
	}
	answer.Alignment = append(answer.Alignment, CalcEntireBlock(align.RSeq, align.QSeq)...)
	return answer
}
