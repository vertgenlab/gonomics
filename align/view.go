package align

import (
	"bytes"
	"fmt"
	"github.com/vertgenlab/gonomics/dna"
	"log"
)

func colTypeToRune(a ColType) rune {
	switch a {
	case ColM:
		return 'M'
	case ColI:
		return 'I'
	case ColD:
		return 'D'
	default:
		log.Fatalf(fmt.Sprintf("Error: unexpected value when converting colType to rune %d", a))
		return '?'
	}
}

func PrintCigar(operations []Cigar) string {
	var buffer bytes.Buffer
	for _, curr := range operations {
		buffer.WriteString(fmt.Sprintf("%d", curr.RunLength))
		buffer.WriteRune(colTypeToRune(curr.Op))
	}
	return buffer.String()
}

func View(alpha []dna.Base, beta []dna.Base, operations []Cigar) string {
	var seqOne, seqTwo bytes.Buffer
	var i, j int
	var count int64
	for _, operation := range operations {
		for count = 0; count < operation.RunLength; count++ {
			switch operation.Op {
			case ColM:
				seqOne.WriteRune(dna.BaseToRune(alpha[i]))
				seqTwo.WriteRune(dna.BaseToRune(beta[j]))
				i, j = i+1, j+1
			case ColI:
				seqOne.WriteRune('-')
				seqTwo.WriteRune(dna.BaseToRune(beta[j]))
				j++
			case ColD:
				seqOne.WriteRune(dna.BaseToRune(alpha[i]))
				seqTwo.WriteRune('-')
				i++
			}
		}
	}
	return seqOne.String() + "\n" + seqTwo.String() + "\n"
}

func LocalView(alpha []dna.Base, beta []dna.Base, operations []Cigar, maxI int64) string {
	var seqOne, seqTwo bytes.Buffer
	var i, j int
	var count int64
	var alignLen int64
	for k := 0; k < len(operations); k++ {
		alignLen += operations[k].RunLength
	}
	var answer []Cigar
	var startCig int64
	var endCig int64
	endCig = int64(len(alpha)) - int64(maxI)
	startCig = int64(len(alpha)) - alignLen - endCig
	fmt.Println("start", startCig, "end", endCig)
	answer = append(answer, Cigar{RunLength: startCig, Op: ColD})
	for x := 0; x < len(operations); x++ {
		answer = append(answer, operations[x])
	}
	answer = append(answer, Cigar{RunLength: endCig, Op: ColD})
	if int(alignLen) != len(alpha) {
		for _, y := range answer {
			for count = 0; count < y.RunLength; count++ {
				switch y.Op {
				case ColM:
					seqOne.WriteRune(dna.BaseToRune(alpha[i]))
					seqTwo.WriteRune(dna.BaseToRune(beta[j]))
					i, j = i+1, j+1
				case ColI:
					seqOne.WriteRune('-')
					seqTwo.WriteRune(dna.BaseToRune(beta[j]))
					j++
				case ColD:
					seqOne.WriteRune(dna.BaseToRune(alpha[i]))
					seqTwo.WriteRune('-')
					i++
				}
			}
		}
	} else {
		for _, operation := range operations {
			for count = 0; count < operation.RunLength; count++ {
				switch operation.Op {
				case ColM:
					seqOne.WriteRune(dna.BaseToRune(alpha[i]))
					seqTwo.WriteRune(dna.BaseToRune(beta[j]))
					i, j = i+1, j+1
				case ColI:
					seqOne.WriteRune('-')
					seqTwo.WriteRune(dna.BaseToRune(beta[j]))
					j++
				case ColD:
					seqOne.WriteRune(dna.BaseToRune(alpha[i]))
					seqTwo.WriteRune('-')
					i++
				}
			}
		}

	}
	return seqOne.String() + "\n" + seqTwo.String() + "\n"
}