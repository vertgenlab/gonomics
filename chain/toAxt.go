package chain

import (
	"github.com/vertgenlab/gonomics/axt"
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
		if each.TDiff > 0 {
			answer.RSeq = append(answer.RSeq, getSequence(target, tIndex, each.TDiff)...)
			answer.QSeq = append(answer.QSeq, dna.CreateAllGaps(int64(each.TDiff))...)
			tIndex += each.TDiff
		}
		if each.QDiff > 0 {
			answer.QSeq = append(answer.QSeq, getSequence(query, qIndex, each.QDiff)...)
			answer.RSeq = append(answer.RSeq, dna.CreateAllGaps(int64(each.QDiff))...)
			qIndex += each.QDiff
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
