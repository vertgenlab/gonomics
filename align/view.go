package align

import (
	"bytes"
	"fmt"

	"github.com/vertgenlab/gonomics/cigar"
	"github.com/vertgenlab/gonomics/dna"
)

// PrintCigar returns the slice of cigar operations as a human-readable string
func PrintCigar(operations []cigar.Cigar) string {
	var buffer bytes.Buffer
	for _, curr := range operations {
		buffer.WriteString(fmt.Sprintf("%d", curr.RunLength))
		buffer.WriteByte(curr.Op)
	}
	return buffer.String()
}

// View takes two sequences and a cigar describing their alignment and returns a
// human-readable alignment of the two sequences.
func View(alpha []dna.Base, beta []dna.Base, operations []cigar.Cigar) string {
	var seqOne, seqTwo bytes.Buffer
	var i, j int
	var count int
	for _, operation := range operations {
		for count = 0; count < operation.RunLength; count++ {
			switch operation.Op {
			case cigar.Match:
				seqOne.WriteRune(dna.BaseToRune(alpha[i]))
				seqTwo.WriteRune(dna.BaseToRune(beta[j]))
				i, j = i+1, j+1
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
	return seqOne.String() + "\n" + seqTwo.String() + "\n"
}

// LocalView returns a human-readable local alignment of two DNA sequences (alpha, beta)
// given the cigar of that alignment and the last aligning base position in alpha.
func LocalView(alpha []dna.Base, beta []dna.Base, operations []cigar.Cigar, maxI int) string {
	var seqOne, seqTwo bytes.Buffer
	var i, j int
	var count int
	var alignLen int
	for k := 0; k < len(operations); k++ {
		alignLen += operations[k].RunLength
	}
	var startCig int
	var endCig int
	endCig = len(alpha) - maxI
	startCig = len(alpha) - alignLen - endCig
	if startCig != 0 {
		operations = append([]cigar.Cigar{{RunLength: startCig, Op: cigar.Deletion}}, operations...)
	}
	if endCig != 0 {
		operations = append(operations, cigar.Cigar{RunLength: endCig, Op: cigar.Deletion})
	}
	fmt.Println("start", startCig, "end", endCig)
	for _, operation := range operations {
		for count = 0; count < operation.RunLength; count++ {
			switch operation.Op {
			case cigar.Match:
				seqOne.WriteRune(dna.BaseToRune(alpha[i]))
				seqTwo.WriteRune(dna.BaseToRune(beta[j]))
				i, j = i+1, j+1
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
	return seqOne.String() + "\n" + seqTwo.String() + "\n"
}
