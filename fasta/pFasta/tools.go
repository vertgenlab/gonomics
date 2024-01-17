package pFasta

import (
	"log"

	"golang.org/x/exp/rand"

	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/dna/pDna"
	"github.com/vertgenlab/gonomics/fasta"
)

// Extract returns a new pFa that is a subsequence of the input pFa, defined by a 
// start (inclusive) and end (exclusive) position, like in bed; makes memory copy
func Extract(input PFasta, start int, end int) PFasta {
	if start >= end {
		log.Fatalf("Error: start must be less than end\n")
	} else if start < 0 || end > len(input.Seq) {
		log.Fatalf("Error: positions out of range\n")
	}
	
	var answer = PFasta{Name: input.Name, Seq: make([]pDna.Float32Base, end-start)}

	for inputIdx := start; inputIdx < end; inputIdx++ {
		answer.Seq[inputIdx-start] = input.Seq[inputIdx]
	}

	return answer
}

// ExtractBed returns a new pFa that is a subsequence of the input pFa
// defined by the bed region
func ExtractBed(input []PFasta, region bed.Bed) PFasta{
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
	} else if len(input[regionIdx].Seq) < region.ChromEnd-region.ChromStart {
		log.Fatalf("Error: region out of range\n")
	}

	var answer = PFasta{Name: input[regionIdx].Name, Seq: make([]pDna.Float32Base, region.ChromEnd-region.ChromStart)}

	for inputIdx := region.ChromStart; inputIdx < region.ChromEnd; inputIdx++ {
		answer.Seq[inputIdx-region.ChromStart] = input[regionIdx].Seq[inputIdx]
	}

	return answer
}

// Sample returns a new Fasta sampled from the given pFasta probability distribution
func Sample(input PFasta) fasta.Fasta {
	var answer = fasta.Fasta{Name: input.Name, Seq: make([]dna.Base, len(input.Seq))}

	var currRand float32 // rand number between 0-1
	for inputIdx := range input.Seq {
		currRand = rand.Float32()
		if currRand < input.Seq[inputIdx].A {
			answer.Seq[inputIdx] = dna.A
		} else if currRand < (input.Seq[inputIdx].C+input.Seq[inputIdx].A) {
			answer.Seq[inputIdx] = dna.C
		} else if currRand < (input.Seq[inputIdx].G+input.Seq[inputIdx].C+input.Seq[inputIdx].A) {
			answer.Seq[inputIdx] = dna.G
		} else {
			answer.Seq[inputIdx] = dna.T
		}
	}

	return answer
}
