package chain

import (
	"log"

	"github.com/vertgenlab/gonomics/axt"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/numbers"
)

// AllToAxt converts a slice of chain structs into a slice of axt structs, taking in the target and query faMaps.
func AllToAxt(ch []Chain, target map[string][]dna.Base, query map[string][]dna.Base) []axt.Axt {
	var answer = make([]axt.Axt, 0)
	for i := range ch {
		answer = append(answer, ToAxt(ch[i], target[ch[i].TName], query[ch[i].QName]))
	}
	return answer
}

// ToAxt converts a chain struct to an axt struct, taking in the target and query chromosome sequences.
func ToAxt(ch Chain, target []dna.Base, query []dna.Base) axt.Axt {
	var tLen, qLen int = ch.TEnd - ch.TStart, ch.QEnd - ch.QStart

	var answer axt.Axt = axt.Axt{
		RName:      ch.TName,
		RStart:     ch.TStart + 1,
		REnd:       ch.TEnd,
		QName:      ch.QName,
		QStart:     ch.QStart + 1,
		QEnd:       ch.QEnd,
		QStrandPos: ch.QStrand,
		Score:      ch.Score,
		RSeq:       make([]dna.Base, 0, tLen),
		QSeq:       make([]dna.Base, 0, qLen),
	}

	var targetFa []dna.Base = make([]dna.Base, len(target))
	copy(targetFa, target)

	var queryFa []dna.Base = make([]dna.Base, len(query))
	copy(queryFa, query)

	switch true {
	case !ch.TStrand && ch.QStrand:
		//axt aligns to target so if target is minus, you need to reverse comp both
		dna.ReverseComplement(targetFa)
		dna.ReverseComplement(queryFa)
	case ch.TStrand && !ch.QStrand:
		//simple case where only query needs to be reverse comp
		dna.ReverseComplement(queryFa)
	case !ch.TStrand && !ch.QStrand:
		//both are negative, you have to swap the target, meaning both get swap, but we are getting the sequence from the fasta you swap query one more time, so you only need to handle target
		dna.ReverseComplement(targetFa)

	default:
		//if they are both positive we do not run the switch
	}

	var tIndex, qIndex int = ch.TStart, ch.QStart
	for _, each := range ch.Alignment {
		answer.RSeq = append(answer.RSeq, getSequence(targetFa, tIndex, each.Size)...)
		answer.QSeq = append(answer.QSeq, getSequence(queryFa, qIndex, each.Size)...)

		tIndex, qIndex = tIndex+each.Size, qIndex+each.Size
		//update idx
		if each.TBases > 0 {
			answer.RSeq = append(answer.RSeq, getSequence(targetFa, tIndex, each.TBases)...)
			answer.QSeq = append(answer.QSeq, dna.CreateAllGaps(each.TBases)...)
			tIndex += each.TBases
		}
		if each.QBases > 0 {
			answer.QSeq = append(answer.QSeq, getSequence(queryFa, qIndex, each.QBases)...)
			answer.RSeq = append(answer.RSeq, dna.CreateAllGaps(each.QBases)...)
			qIndex += each.QBases
		}
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
	var answer int = 0
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
	var target, query int = 0, 0
	for i := 0; i < len(rSeq) && (rSeq[i] == dna.Gap || qSeq[i] == dna.Gap); i++ {
		if rSeq[i] == dna.Gap {
			query++
		}
		if qSeq[i] == dna.Gap {
			target++
		}
	}
	return target, query
}

func CalcEntireBlock(rSeq []dna.Base, qSeq []dna.Base) []BaseStats {
	var answer []BaseStats
	var curr BaseStats
	if len(rSeq) != len(qSeq) {
		log.Fatalf("Error input sequences should match in length\n")
	}
	for i := 0; i < len(rSeq) && i < len(qSeq); {
		//first count matching seq
		curr = BaseStats{
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
			i += numbers.Max(curr.TBases, curr.QBases)
		}
	}
	return answer
}

func AxtToChain(align *axt.Axt, tLen int, qLen int, id int) Chain {
	var answer Chain = Chain{
		Score:     int(align.Score),
		TName:     align.RName,
		TSize:     tLen,
		TStrand:   true,
		TStart:    int(align.RStart) - 1,
		TEnd:      int(align.REnd),
		QName:     align.QName,
		QSize:     qLen,
		QStrand:   align.QStrandPos,
		QStart:    int(align.QStart) - 1,
		QEnd:      int(align.QEnd),
		Alignment: make([]BaseStats, 0),
		Id:        id,
	}
	answer.Alignment = append(answer.Alignment, CalcEntireBlock(align.RSeq, align.QSeq)...)
	return answer
}
