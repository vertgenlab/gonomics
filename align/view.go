package align

import (
	"bytes"
	"fmt"
	"github.com/craiglowe/gonomics/common"
	"github.com/craiglowe/gonomics/dna"
)

func colTypeToRune(a colType) rune {
	switch a {
	case colM:
		return 'M'
	case colI:
		return 'I'
	case colD:
		return 'D'
	default:
		common.Exit(fmt.Sprintf("Error: unexpected value when converting colType to rune %d", a))
		return '?'
	}
}

func printCigar(operations []cigar) string {
	var buffer bytes.Buffer
	for _, curr := range operations {
		buffer.WriteString(fmt.Sprintf("%d", curr.runLength))
		buffer.WriteRune(colTypeToRune(curr.op))
	}
	return buffer.String()
}

func View(alpha []dna.Base, beta []dna.Base, operations []cigar) string {
	var seqOne, seqTwo bytes.Buffer
	var i, j int
	var count int64
	for _, operation := range operations {
		for count = 0; count < operation.runLength; count++ {
			switch operation.op {
			case colM:
				seqOne.WriteRune(dna.BaseToRune(alpha[i]))
				seqTwo.WriteRune(dna.BaseToRune(beta[j]))
				i, j = i+1, j+1
			case colI:
				seqOne.WriteRune('-')
				seqTwo.WriteRune(dna.BaseToRune(beta[j]))
				j++
			case colD:
				seqOne.WriteRune(dna.BaseToRune(alpha[i]))
				seqTwo.WriteRune('-')
				i++
			}
		}
	}
	return seqOne.String() + "\n" + seqTwo.String() + "\n"
}
