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
func Extract(input []PFasta, start int, end int, outputName string, chrom string) PFasta {

	//else if input.Name != chrom {
	//	log.Fatalf("Error: input sequence name does not match requested chrom.")
	//}

	regionInInput := false
	regionIdx := 0

	for inputIdx, inputpFa := range input {
		if inputpFa.Name == chrom {
			regionInInput = true
			regionIdx = inputIdx
			break
		}
	}

	if start >= end {
		log.Fatalf("Error: start must be less than end\n")
	} else if start < 0 || end > len(input[regionIdx].Seq) {
		log.Fatalf("Error: positions out of range\n")
	}

	if !regionInInput {
		log.Fatalf("Error: region not in input\n")
	}

	var outName string
	if len(outputName) > 0 {
		outName = outputName
	} else {
		outName = chrom
	}

	var answer = PFasta{Name: outName, Seq: make([]pDna.Float32Base, end-start)}

	for inputIdx := start; inputIdx < end; inputIdx++ {
		answer.Seq[inputIdx-start] = input[regionIdx].Seq[inputIdx]
	}

	return answer
}

// ExtractBed returns a new pFa that is a subsequence of the input pFa
// defined by the bed region
func ExtractBed(input []PFasta, region bed.Bed, outputName string) PFasta {
	return Extract(input, region.ChromStart, region.ChromEnd, outputName, region.Chrom)
}

// Sample returns a new Fasta sampled from the given pFasta probability distribution

// //// ASK THIS QUESTION TOMORROW: do we want to be able to specify which sequence in the input we want a sample of?
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
