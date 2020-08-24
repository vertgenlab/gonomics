package giraf

import (
	"github.com/vertgenlab/gonomics/cigar"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/simpleGraph"
	"log"
)

func GirafToExplicitCigar(giraf *Giraf, graph *simpleGraph.SimpleGraph) []*cigar.Cigar {
	var answer []*cigar.Cigar
	var seqIdx, refIdx, pathIdx int
	refIdx = giraf.Path.TStart
	var k, runLenCount int64

	if giraf.Aln[0].Op == '*' {
		return nil
	}

	for i := 0; i < len(giraf.Aln); i++ {
		switch giraf.Aln[i].Op {
		case 'M':
			runLenCount = 0

			for k = 0; k < giraf.Aln[i].RunLength; k++ {
				if refIdx > len(graph.Nodes[giraf.Path.Nodes[pathIdx]].Seq)-1 {
					pathIdx++
					refIdx = 0
				}
				if giraf.Seq[seqIdx] == graph.Nodes[giraf.Path.Nodes[pathIdx]].Seq[refIdx] {
					runLenCount++
				} else {
					if runLenCount > 0 {
						// Append the matching bases so far
						answer = append(answer, &cigar.Cigar{RunLength: runLenCount, Op: '='})
					}
					// Append the mismatch base
					if answer == nil {
						answer = append(answer, &cigar.Cigar{RunLength: 1, Op: 'X', Sequence: []dna.Base{giraf.Seq[k]}})
					} else if answer[len(answer)-1].Op == 'X' {
						answer[len(answer)-1].RunLength++
						answer[len(answer)-1].Sequence = append(answer[len(answer)-1].Sequence, giraf.Seq[k])
					} else {
						answer = append(answer, &cigar.Cigar{RunLength: 1, Op: 'X', Sequence: []dna.Base{giraf.Seq[k]}})
					}
					runLenCount = 0
				}
				seqIdx++
				refIdx++
			}

			if runLenCount > 0 {
				answer = append(answer, &cigar.Cigar{RunLength: runLenCount, Op: '='})
			}

		case 'I':
			var insSeq []dna.Base
			for k = 0; k < giraf.Aln[i].RunLength; k++ {
				insSeq = append(insSeq, giraf.Seq[seqIdx])
				seqIdx++
			}
			answer = append(answer, &cigar.Cigar{RunLength: giraf.Aln[i].RunLength, Op: 'I', Sequence: insSeq})

		case 'X':
			log.Println("WARNING: The input cigar already has explicit formatting")
			return giraf.Aln

		case '=':
			log.Println("WARNING: The input cigar already has explicit formatting")
			return giraf.Aln

		default:
			answer = append(answer, giraf.Aln[i])
			if cigar.ConsumesReference(giraf.Aln[i].Op) {
				refIdx += int(giraf.Aln[i].RunLength)
			}
			if cigar.ConsumesQuery(giraf.Aln[i].Op) {
				seqIdx += int(giraf.Aln[i].RunLength)
			}
		}
	}
	return answer
}
