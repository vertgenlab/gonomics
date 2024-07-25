package genomeGraph

import (
	"strings"

	"github.com/vertgenlab/gonomics/cigar"
	"github.com/vertgenlab/gonomics/dna"
)

func View(alpha []dna.Base, beta []dna.Base, operations []cigar.Cigar, startI, endI, startJ, endJ int) string {
	var seqOne, seqTwo strings.Builder
	var i, j = 0, 0
	// Print leading bases in alpha
	for ; i < startI; i++ {
		seqOne.WriteRune(dna.BaseToRune(alpha[i]))
		seqTwo.WriteRune('-')
	}

	// Print leading bases in beta
	for ; j < startJ; j++ {
		seqOne.WriteRune('-')
		seqTwo.WriteRune(dna.BaseToRune(beta[j]))
	}

	// Print the aligned region
	for _, operation := range operations {
		for count := 0; count < operation.RunLength && i <= endI && j <= endJ; count++ {
			switch operation.Op {
			case cigar.Match, cigar.Equal, cigar.Mismatch:
				seqOne.WriteRune(dna.BaseToRune(alpha[i]))
				seqTwo.WriteRune(dna.BaseToRune(beta[j]))
				i++
				j++
			case cigar.Insertion:
				seqOne.WriteRune('-')
				seqTwo.WriteRune(dna.BaseToRune(beta[j]))
				j++
			case cigar.Deletion:
				seqOne.WriteRune(dna.BaseToRune(alpha[i]))
				seqTwo.WriteRune('-')
				i++
			}
		}
	}

	// Print trailing bases in alpha (if any)
	for ; i < endI; i++ {
		seqOne.WriteRune(dna.BaseToRune(alpha[i]))
		seqTwo.WriteRune('-')
	}
	seqOne.WriteByte('\n')

	// Print trailing bases in beta (if any)
	for ; j < endJ; j++ {
		seqOne.WriteRune('-')
		seqTwo.WriteRune(dna.BaseToRune(beta[j]))
	}
	seqTwo.WriteByte('\n')

	return seqOne.String() + seqTwo.String()
}
