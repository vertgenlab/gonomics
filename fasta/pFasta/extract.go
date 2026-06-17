package pFasta

import (
	"fmt"
	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/dna/pDna"
	"log"
)

// Extract returns a new pFa that is a subsequence of the input pFa, defined by a
// start (inclusive) and end (exclusive) position, like in bed; makes memory copy
func Extract(input []PFasta, start int, end int, outputName string, chrom string, takeCoords bool) PFasta {

	chromIdx := checkIfChromInPfasta(input, chrom)

	if start >= end {
		log.Fatalf("Error: start must be less than end\n")
	} else if start < 0 || end > len(input[chromIdx].Seq) {
		log.Fatalf("Error: positions out of range\n")
	}

	var outName string
	if takeCoords {
		outName = fmt.Sprintf("%s:%v-%v", chrom, start, end)
	} else if len(outputName) > 0 {
		outName = outputName
	} else {
		outName = chrom
	}

	var answer = PFasta{Name: outName, Seq: make([]pDna.Float32Base, end-start)}

	for inputIdx := start; inputIdx < end; inputIdx++ {
		answer.Seq[inputIdx-start] = input[chromIdx].Seq[inputIdx]
	}

	return answer
}

// ExtractBed returns a pFa that has a list of subsequences of the input pFa
// defined by the regions in the bed region
// takeCoords specifies if name fields in output should be original names in region or identified by ChromStart and ChromEnd
func ExtractBed(input []PFasta, region []bed.Bed, takeCoords bool) []PFasta {
	answer := make([]PFasta, 0)
	for _, reg := range region {
		answer = append(answer, Extract(input, reg.ChromStart, reg.ChromEnd, "", reg.Chrom, takeCoords))
	}
	return answer
}
