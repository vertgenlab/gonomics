package cigar

import (
	"fmt"
	"github.com/vertgenlab/gonomics/dna"
	"log"
)

//TODO: build functions to rip linear sequence for cigar from []*Fasta or a *SimpleGraph

// Make an existing sam-formatted cigar include the explicit sequence so that the read sequence is no longer necessary
// Mismatch character 'X'
// Sequence additions are lowercase in an effort to improve human readability
func MakeExplicit(cigar []*Cigar, sequence []dna.Base, reference []dna.Base) []*Cigar {
	fmt.Println("Started")
	var answer []*Cigar
	var seqIdx, refIdx int
	var k, runLenCount int64

	if cigar[0].Op == '*' {
		return []*Cigar{&Cigar{Op: '*'}}
	}

	for i := 0; i < len(cigar); i++ {
		switch cigar[i].Op {
		case 'M':
			runLenCount = 0
			for k = 0; k < cigar[i].RunLength; k++ {
				if sequence[seqIdx] == reference[refIdx] {
					runLenCount++
				} else {
					if runLenCount > 0 {
						// Append the matching bases so far
						answer = append(answer, &Cigar{RunLength: runLenCount, Op: '='})
					}
					// Append the mismatch base
					if answer == nil {
						answer = append(answer, &Cigar{RunLength: 1, Op: 'X', Sequence: []dna.Base{sequence[k]}})
					} else if answer[len(answer) - 1].Op == 'X' {
						answer[len(answer) - 1].RunLength++
						answer[len(answer) - 1].Sequence = append(answer[len(answer) - 1].Sequence ,sequence[k])
					} else {
						answer = append(answer, &Cigar{RunLength: 1, Op: 'X', Sequence: []dna.Base{sequence[k]}})
					}
					runLenCount = 0
				}
				seqIdx++
				refIdx++
			}

			if runLenCount > 0 {
				answer = append(answer, &Cigar{RunLength: runLenCount, Op: '='})
			}

		case 'I':
			var insSeq []dna.Base
			for k = 0; k < cigar[i].RunLength; k++ {
				insSeq = append(insSeq, sequence[seqIdx])
				seqIdx++
			}
			answer = append(answer, &Cigar{RunLength: cigar[i].RunLength, Op: 'I', Sequence: insSeq})

		case 'X':
			log.Println("WARNING: The input cigar already has explicit formatting")
			return cigar

		case '=':
			log.Println("WARNING: The input cigar already has explicit formatting")
			return cigar

		default:
			answer = append(answer, cigar[i])
			if ConsumesReference(cigar[i].Op) {
				refIdx += int(cigar[i].RunLength)
			}
			if ConsumesQuery(cigar[i].Op) {
				seqIdx += int(cigar[i].RunLength)
			}
		}
	}
	return answer
}

// Imports an explicit cigar and a reference sequence (startpos of ref = startpos of cigar)
// and constructs the entire read sequence from the cigar
func GetExplicitSequence(cigar []*Cigar, reference []dna.Base) []dna.Base {
	var answer []dna.Base
	var refIdx int

	for i := 0; i < len(cigar); i++ {
		switch cigar[i].Op {
		case 'I':
			answer = append(answer, cigar[i].Sequence...)
		case 'D':
			refIdx += int(cigar[i].RunLength)
		case 'X':
			answer = append(answer, cigar[i].Sequence...)
			refIdx += int(cigar[i].RunLength)
		default:
			answer = append(answer, reference[refIdx:refIdx+int(cigar[i].RunLength)]...)
			if ConsumesReference(cigar[i].Op) {
				refIdx += int(cigar[i].RunLength)
			}
		}
	}
	return answer
}
