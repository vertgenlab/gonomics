package genomeGraph

import (
	"bytes"
	"github.com/vertgenlab/gonomics/cigar"
	"github.com/vertgenlab/gonomics/common"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/sam"
	"strings"
)

func calcExtension(seq []dna.Base) int64 {
	var maxScore int64 = 0
	for bases := 0; bases < len(seq); bases++ {
		maxScore += HumanChimpTwoScoreMatrix[seq[bases]][seq[bases]]
	}
	return maxScore/600 + int64(len(seq))
}

func LocalView(samLine *sam.Aln, ref []*Node) string {
	var seqOne, seqTwo bytes.Buffer

	var operations []*cigar.Cigar = samLine.Cigar
	var i int64 = int64(samLine.Pos) - 1
	var j int64 = 0
	var count int
	words := strings.Split(samLine.RName, "_")
	var alpha []dna.Base = ref[common.StringToInt64(words[1])].Seq
	var beta []dna.Base = samLine.Seq

	for _, operation := range operations {
		for count = 0; count < operation.RunLength; count++ {
			switch operation.Op {
			case 'M':
				seqOne.WriteRune(dna.BaseToRune(alpha[i]))
				seqTwo.WriteRune(dna.BaseToRune(beta[j]))
				i, j = i+1, j+1
			case 'I':
				seqOne.WriteRune('-')
				seqTwo.WriteRune(dna.BaseToRune(beta[j]))
				j++
			case 'D':
				seqOne.WriteRune(dna.BaseToRune(alpha[i]))
				seqTwo.WriteRune('-')
				i++
			case 'S':
				seqOne.WriteRune('-')
				seqTwo.WriteRune('-')
			}
		}
	}
	return cigar.ToString(operations) + "\n" + seqOne.String() + "\n" + seqTwo.String() + "\n"
}
