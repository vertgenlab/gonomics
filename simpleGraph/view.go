package simpleGraph

import (
	"bytes"
	"github.com/vertgenlab/gonomics/cigar"
	"github.com/vertgenlab/gonomics/common"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/sam"
	"github.com/vertgenlab/gonomics/giraf"
	"strings"
	"fmt"
)

func calcExtension(seq []dna.Base) int64 {
	var maxScore int64 = 0
	for bases := 0; bases < len(seq); bases++ {
		maxScore += HumanChimpTwoScoreMatrix[seq[bases]][seq[bases]]
	}
	return maxScore/600 + int64(len(seq))
}

func LocalView(samLine *sam.SamAln, ref []*Node) string {
	var seqOne, seqTwo bytes.Buffer
	var operations []*cigar.Cigar = samLine.Cigar
	var i int64 = samLine.Pos - 1
	var j int64 = 0
	var count int64
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


func GirafPathToSeq(p *giraf.Path, gg *SimpleGraph) []dna.Base {
	var answer []dna.Base
	if len(p.Nodes) == 0 {
		return nil
	} else if len(p.Nodes) == 1 {
		answer = append( answer, gg.Nodes[p.Nodes[0]].Seq[p.TStart:p.TEnd]...)
	} else if len(p.Nodes) == 2 {
		answer = append(append(answer, gg.Nodes[p.Nodes[0]].Seq[p.TStart:]...), gg.Nodes[p.Nodes[1]].Seq[:p.TEnd]...)
	} else {
		answer = append(answer, gg.Nodes[p.Nodes[0]].Seq[p.TStart:]...)
		for i := 0; i < len(p.Nodes)-1; i++ {
			answer = append(answer, gg.Nodes[p.Nodes[i]].Seq...)
		}
		answer = append(answer, gg.Nodes[p.Nodes[len(p.Nodes)-1]].Seq...)
	}
	return answer
}

func ViewGirafAlign(alignment *giraf.Giraf, genome *SimpleGraph) string {
	if alignment.Path == nil {
		return fmt.Sprintf("Unmapped Alignment:\n%s\n", giraf.GirafToString(alignment))
	} else {
		var seqOne, seqTwo bytes.Buffer
		var operations []*cigar.Cigar = alignment.Aln
		var i int = alignment.Path.TStart
		var j int = alignment.QStart
		var count int64
		var alpha []dna.Base = GirafPathToSeq(alignment.Path, genome)
		var beta []dna.Base = alignment.Seq
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
		return fmt.Sprintf("%s\n%s\n", seqOne.String(), seqTwo.String())
	}
}
