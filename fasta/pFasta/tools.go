package pFasta

import (
	"log"

	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/dna/pDna"
	"github.com/vertgenlab/gonomics/fasta"
	"math/rand"
)

// Extract returns a new pFa that is a subsequence of the input pFa, defined by a
// start (inclusive) and end (exclusive) position, like in bed; makes memory copy
func Extract(input PFasta, start int, end int, outputName string) PFasta {
	if start >= end {
		log.Fatalf("Error: start must be less than end\n")
	} else if start < 0 || end > len(input.Seq) {
		log.Fatalf("Error: positions out of range\n")
	}

	var outName string
	if len(outputName) > 0 {
		outName = outputName
	} else {
		outName = input.Name
	}

	var answer = PFasta{Name: outName, Seq: make([]pDna.Float32Base, end-start)}

	for inputIdx := start; inputIdx < end; inputIdx++ {
		answer.Seq[inputIdx-start] = input.Seq[inputIdx]
	}

	return answer
}

// ExtractBed returns a new pFa that is a subsequence of the input pFa
// defined by the bed region
func ExtractBed(input []PFasta, region bed.Bed, outputName string) PFasta {
	regionInInput := false
	regionIdx := 0
	for inputIdx, inputpFa := range input {
		if inputpFa.Name == region.Chrom {
			regionInInput = true
			regionIdx = inputIdx
			break
		}
	}

	if !regionInInput {
		log.Fatalf("Error: region not in input\n")
	}

	return Extract(input[regionIdx], region.ChromStart, region.ChromEnd, outputName)
}

// Sample returns a new Fasta sampled from the given pFasta probability distribution
func Sample(input PFasta) fasta.Fasta {
	var answer = fasta.Fasta{Name: input.Name, Seq: make([]dna.Base, len(input.Seq))}
	var currRand float32
	for inputIdx := range input.Seq {
		currRand = rand.Float32()
		if currRand < input.Seq[inputIdx].A {
			answer.Seq[inputIdx] = dna.A
		} else if currRand < (input.Seq[inputIdx].C + input.Seq[inputIdx].A) {
			answer.Seq[inputIdx] = dna.C
		} else if currRand < (input.Seq[inputIdx].G + input.Seq[inputIdx].C + input.Seq[inputIdx].A) {
			answer.Seq[inputIdx] = dna.G
		} else {
			answer.Seq[inputIdx] = dna.T
		}
	}

	return answer
}
