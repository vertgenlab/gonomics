package qDna

import (
	"github.com/vertgenlab/gonomics/align"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/fastq"
	"github.com/vertgenlab/gonomics/sam"
	"reflect"
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

func UngappedAlignLen(cig []align.Cigar) int64 {
	var reds int64
	for i := 0; i < len(cig); i++ {
		if cig[i].Op != align.ColD {
			reds = reds + cig[i].RunLength
		}
	}
	return reds
}

func UngappedQueryLen(cig []align.Cigar) int64 {
	var reds int64
	for i := 0; i < len(cig); i++ {
		if cig[i].Op != align.ColI {
			reds = reds + cig[i].RunLength
		}
	}
	return reds
}

func GSW(ref []*QFrag, reads []*fastq.Fastq) []*sam.SamAln {
	var answer []*sam.SamAln = make([]*sam.SamAln, len(reads))
	var reverseRead []*QBase
	//var query []*qDna.QBase
	//var reverse []*qDna.QBase
	var score float64
	var alignment []align.Cigar
	var lowRef, lowQuery, highQuery int64
	var reverseFastq *fastq.Fastq

	var bestScore float64 = 0
	//var minI, minJ, maxJ int64
	//var qualBase []rune
	//var sequence string
	//var flag int64
	var i, j int
	var currRead []*QBase
	for i = 0; i < len(reads); i++ {
		currRead = FromFastq(reads[i])
		reverseFastq = fastq.ReverseComplementFastq(reads[i])
		reverseRead = FromFastq(reverseFastq)

		var currBest sam.SamAln = sam.SamAln{QName: reads[i].Name, Flag: 0, RName: "", Pos: 0, MapQ: 255, Cigar: "", RNext: "*", PNext: 0, TLen: 0, Seq: "", Qual: "", Extra: ""}
		bestScore = 0
		for j = 0; j < len(ref); j++ {

			//query := FromDna(reads[j].Seq, ErrorRate(reads[j].Qual))

			score, alignment, lowRef, _, lowQuery, highQuery = SmithWaterman(ref[j].Seq, currRead, HumanChimpTwoScoreMatrix, -600)
			if score > bestScore {
				bestScore = score

				currBest.Flag = 0
				currBest.RName = ref[j].From[0].Chr
				currBest.Pos = lowRef
				currBest.MapQ = 255
				currBest.Cigar = align.PrintCigar(alignment)
				currBest.Seq = dna.BasesToString(reads[i].Seq[lowQuery:highQuery])
				currBest.Qual = string(reads[i].Qual[lowQuery:highQuery])

			}
			score, alignment, lowRef, _, lowQuery, highQuery = SmithWaterman(ref[j].Seq, reverseRead, HumanChimpTwoScoreMatrix, -600)
			if score > bestScore {
				bestScore = score

				currBest.Flag = 16
				currBest.RName = ref[j].From[0].Chr
				currBest.Pos = lowRef
				currBest.MapQ = 255
				currBest.Cigar = align.PrintCigar(alignment)
				currBest.Seq = dna.BasesToString(reverseFastq.Seq[lowQuery:highQuery])
				currBest.Qual = string(reverseFastq.Qual[lowQuery:highQuery])

			}
			answer[i] = &currBest
			//answer = append(answer, &sam.SamAln{QName: qName, Flag: flag, RName: readName, Pos: minI, MapQ: mappingQ, Cigar: reds, RNext: rNext, PNext: pNext, TLen: tlen, Seq: sequence, Qual: string(qualBase), Extra: ""})
		}
	}
	return answer
}

/*
func UngappedAlign(alpha []*QBase, beta []*QBase, alphaOffset int, betaOffset int, scoreMatrix [][]float64) float64 {
}*/
