package pFasta

import (
	"fmt"
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

func FaToPfa(fa fasta.Fasta) PFasta {
	var pfa PFasta
	var pbase pDna.Float32Base

	pfa.Name = fa.Name
	for i := range fa.Seq {
		pbase = pDna.DnaToPdna(fa.Seq[i])
		pfa.Seq = append(pfa.Seq, pbase)
	}

	return pfa
}
