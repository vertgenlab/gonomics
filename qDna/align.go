package qDna

import (
	"github.com/vertgenlab/gonomics/align"
	"reflect"
	//"github.com/vertgenlab/gonomics/common"
)

type NoGapAln struct {
	Start int
	End   int
	Score float64
}

//aligner name gsw, graph smith-waterman
// O=600 E=150
var HumanChimpTwoScoreMatrix = [][]float64{
	{90, -330, -236, -356},
	{-330, 100, -318, -236},
	{-236, -318, 100, -330},
	{-356, -236, -330, 90},
}

func QDnaScore(alpha *QBase, beta *QBase, scoreMatrix [][]float64) float64 {
	var sum float64 = 0
	a := reflect.ValueOf(alpha).Elem()
	b := reflect.ValueOf(beta).Elem()
	for x := 0; x < a.NumField(); x++ {
		for y := 0; y < b.NumField(); y++ {
			sum += scoreMatrix[x][y] * a.Field(x).Float() * b.Field(y).Float()
		}
	}
	return sum
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

func reverseCigar(alpha []align.Cigar) {
	for i, j := 0, len(alpha)-1; i < j; i, j = i+1, j-1 {
		alpha[i], alpha[j] = alpha[j], alpha[i]
	}
}

func UngappedAlignLen(cig []align.Cigar) int64{
	var reds int64
	for i := 0; i < len(cig); i++ {
		if cig[i].Op != align.ColD{
			reds = reds + cig[i].RunLength
		}
	}
	return reds
}

func UngappedQueryLen(cig []align.Cigar) int64{
	var reds int64
	for i := 0; i < len(cig); i++ {
		if cig[i].Op != align.ColI{
			reds = reds + cig[i].RunLength
		}
	}
	return reds
}
/*
func UngappedAlign(alpha []*QBase, beta []*QBase, alphaOffset int, betaOffset int, scoreMatrix [][]float64) float64 {
}*/