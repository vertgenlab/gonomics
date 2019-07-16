package fastq

import (
	"github.com/vertgenlab/gonomics/align"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/qDna"
	"github.com/vertgenlab/gonomics/sam"
)

func GSW(ref []*fasta.Fasta, reads []*Fastq) []*sam.SamAln {
	var answer []*sam.SamAln
	var chr *qDna.QFrag
	//var query []*qDna.QBase
	//var reverse []*qDna.QBase
	//var bestScore float64
	var bestCigar []align.Cigar
	var minI, minJ, maxJ int64
	var qualBase []rune
	var sequence string
	var flag int64
	var i, j int64
	for i = int64(0); i < int64(len(ref)); i++ {
		chr = qDna.FromFasta(ref[i])
		for j = int64(0); j < int64(len(reads)); j++ {
			query := FromDna(reads[j].Seq, ErrorRate(reads[j].Qual))
			reverseFastq := ReverseComplementFastq(reads[j])
			reverse := FromFastq(reverseFastq)
			score, alignment, lowI, _, lowJ, highJ := qDna.SmithWaterman(chr.Seq, query, qDna.HumanChimpTwoScoreMatrix, -600)
			reverseScore, negAlignment, rLowI, _, rLowJ, rHighJ := qDna.SmithWaterman(chr.Seq, reverse, qDna.HumanChimpTwoScoreMatrix, -600)
			//sam file fields
			qName := reads[j].Name
			readName := ref[i].Name
			mappingQ := int64(255)
			reds := align.PrintCigar(bestCigar)
			rNext := "*"
			pNext := int64(0)
			tlen := int64(0)
			if score < reverseScore {
				//bestScore = reverseScore
				bestCigar = negAlignment
				minI = rLowI
				//maxI = rHighI
				minJ = rLowJ
				maxJ = rHighJ
				qualBase = reverseFastq.Qual
				sequence = dna.BasesToString(reverseFastq.Seq[minJ:maxJ])
				flag = int64(16)
			} else {
				//bestScore = score
				bestCigar = alignment
				minI = lowI
				//maxI = highI
				minJ = lowJ
				maxJ = highJ
				qualBase = reads[j].Qual
				sequence = dna.BasesToString(reads[j].Seq[minJ:maxJ])
				flag = int64(0)
			}
			answer = append(answer, &sam.SamAln{QName: qName, Flag: flag, RName: readName, Pos: minI, MapQ: mappingQ, Cigar: reds, RNext: rNext, PNext: pNext, TLen: tlen, Seq: sequence, Qual: string(qualBase), Extra: ""})
		}
	}
	return answer
}