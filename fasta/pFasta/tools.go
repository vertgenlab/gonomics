package pFasta

import (
	"log"

	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/dna/pDna"
	"github.com/vertgenlab/gonomics/fasta"
	"math/rand"
)

// checks if input pFasta has a sequence with chrom as name and returns its index
func checkIfChromInPfasta(input []PFasta, chrom string) int {
	chromInInput := false
	var answer int
	
	for inputIdx, inputpFa := range input {
		if inputpFa.Name == chrom {
			chromInInput = true
			answer = inputIdx
		}
	}

	if !chromInInput {
		log.Fatalf("Error: input sequence name does not match requested chrom.")
	}
	
	return answer
}

// Extract returns a new pFa that is a subsequence of the input pFa, defined by a
// start (inclusive) and end (exclusive) position, like in bed; makes memory copy
func Extract(input []PFasta, start int, end int, outputName string, chrom string) PFasta {

	chromIdx := checkIfChromInPfasta(input, chrom)

	if start >= end {
		log.Fatalf("Error: start must be less than end\n")
	} else if start < 0 || end > len(input[chromIdx].Seq) {
		log.Fatalf("Error: positions out of range\n")
	}

	var outName string
	if len(outputName) > 0 {
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

// ExtractBed returns a new pFa that is a subsequence of the input pFa
// defined by the region in the bed region (currently only accepts one region in bed file)
// TODO: add multiple regions in bed
func ExtractBed(input []PFasta, region bed.Bed, outputName string) PFasta {
	return Extract(input, region.ChromStart, region.ChromEnd, outputName, region.Chrom)
}

// Sample returns a new Fasta sampled from the given pFasta probability distribution
func Sample(input []PFasta, chrom string) fasta.Fasta {
	chromIdx := checkIfChromInPfasta(input, chrom)

	var answer = fasta.Fasta{Name: input[chromIdx].Name, Seq: make([]dna.Base, len(input[chromIdx].Seq))}
	var currRand float32
	for inputIdx := range input[chromIdx].Seq {
		currRand = rand.Float32()
		if currRand < input[chromIdx].Seq[inputIdx].A {
			answer.Seq[inputIdx] = dna.A
		} else if currRand < (input[chromIdx].Seq[inputIdx].C + input[chromIdx].Seq[inputIdx].A) {
			answer.Seq[inputIdx] = dna.C
		} else if currRand < (input[chromIdx].Seq[inputIdx].G + input[chromIdx].Seq[inputIdx].C + input[chromIdx].Seq[inputIdx].A) {
			answer.Seq[inputIdx] = dna.G
		} else {
			answer.Seq[inputIdx] = dna.T
		}
	}

	return answer
}
