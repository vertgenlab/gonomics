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
func Extract(input PFasta, start int, end int, outputName string) PFasta {
	if start >= end {
		log.Fatalf("Error: start must be less than end\n")
	} else if start < 0 || end > len(input.Seq) {
		log.Fatalf("Error: positions out of range\n")
	}

	var answer = PFasta{Name: outputName, Seq: make([]pDna.Float32Base, end-start)}

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

	var currRand float32 // rand number between 0-1
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

// PAlnPosToRefPos is AlnPosToRefPos in fasta/multiFa.go but for pFasta
// Consider using pAlnPosToRefPosCounter instead if tracking refStart and alnStart will be beneficial, e.g. when working through entire chromosomes
func PAlnPosToRefPos(record PFasta, AlnPos int) int {
	return PAlnPosToRefPosCounter(record, AlnPos, 0, 0)
}

// PAlnPosToRefPosCounter is AlnPosToRefPosCounter in fasta/multiFa.go but for pFasta
func PAlnPosToRefPosCounter(record PFasta, AlnPos int, refStart int, alnStart int) int {
	return PAlnPosToRefPosCounterSeq(record.Seq, AlnPos, refStart, alnStart)
}

// PAlnPosToRefPosCounterSeq is AlnPosToRefPosCounterSeq in fasta/multiFa.go but for pFasta
func PAlnPosToRefPosCounterSeq(record []pDna.Float32Base, AlnPos int, refStart int, alnStart int) int {
	if alnStart > AlnPos {
		refStart, alnStart = 0, 0 //in case the alnStart was improperly set (greater than the desired position, we reset the counters to 0.
	}
	for t := alnStart; t < AlnPos; t++ {
		if t == len(record) {
			log.Fatalf("Ran out of chromosome.")
		} else if !pDna.IsGap(record[t]) {
			refStart++
		}
	}
	return refStart
}
