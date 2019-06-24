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

func printCigar(operations []Cigar) string {
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
