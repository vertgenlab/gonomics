package qDna

import (
//"fmt"
//"github.com/vertgenlab/gonomics/common"
)

type NoGapAln struct {
	Start int
	End   int
	Score float64
}

func PairwiseAverage(alpha *QFrag, beta *QFrag, start int64, end int64, name string) *QFrag {

	answer := &QFrag{Seq: nil, From: []*Location{&Location{Assembly: "", Chr: name, Start: start, End: end}}, Fwd: nil, Rev: nil}
	//Max or min, haven't decided if it's necessary.
	//Just trying to handle a potential error when sequences are uneven
	//length = common.Min(len(alpha.Seq), len(beta.Seq))
	for i := 0; i < len(alpha.Seq); i++ {
		tmpA := (alpha.Seq[i].A + beta.Seq[i].A) / 2
		tmpC := (alpha.Seq[i].C + beta.Seq[i].C) / 2
		tmpG := (alpha.Seq[i].G + beta.Seq[i].G) / 2
		tmpT := (alpha.Seq[i].T + beta.Seq[i].T) / 2
		answer.Seq = append(answer.Seq, &QBase{A: tmpA, C: tmpC, G: tmpG, T: tmpT})
	}
	return answer
}

/*
func UngappedAlign(alpha []*QBase, beta []*QBase, alphaOffset int, betaOffset int, scoreMatrix [][]float64) float64 {

}*/
