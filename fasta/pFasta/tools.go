package pFasta

import (
	"log"
	"math/rand"

	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/dna/pDna"
	"github.com/vertgenlab/gonomics/fasta"
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

// FaToPfa returns a pFasta representation of the given Fasta sequence, start inclusive, end exclusive
func FaToPfa(input fasta.Fasta, start int, end int) PFasta {
	if end == -1 {
		end = len(input.Seq)
	} else if end > len(input.Seq) {
		log.Fatalf("Requested end argument (%v) out of range.", end)
	}

	answer := PFasta{Name: input.Name, Seq: make([]pDna.Float32Base, end-start)}

	fasta.ToUpper(input)
	var base dna.Base
	var idx int
	for idx, base = range input.Seq[start:end] {
		if base == dna.A {
			answer.Seq[idx] = pDna.Float32Base{A: 1, C: 0, G: 0, T: 0}
		} else if base == dna.C {
			answer.Seq[idx] = pDna.Float32Base{A: 0, C: 1, G: 0, T: 0}
		} else if base == dna.G {
			answer.Seq[idx] = pDna.Float32Base{A: 0, C: 0, G: 1, T: 0}
		} else if base == dna.T {
			answer.Seq[idx] = pDna.Float32Base{A: 0, C: 0, G: 0, T: 1}
		} else if base == dna.N {
			//answer.Seq[idx] = pDna.Float32Base{A: 0, C: 0, G: 0, T: 0} // Raven's change: I think N should be all 0.25, not all 0
			answer.Seq[idx] = pDna.Float32Base{A: 0.25, C: 0.25, G: 0.25, T: 0.25}
		} else if base == dna.Gap {
			//log.Fatalf("Must specify a sequence without gaps.") // Raven's change: I think gap should be all 0, not illegal
			answer.Seq[idx] = pDna.Float32Base{A: 0, C: 0, G: 0, T: 0}
		}
	}

	return answer
}

// MultiFaToPfa returns a pFasta representation of the given Fasta sequence
func MultiFaToPfa(inputFaFilename string, start int, end int, chrom string) PFasta {
	inputFa := fasta.Read(inputFaFilename)
	chromInInput := false
	var answer PFasta

	if len(inputFa) == 1 {
		if chrom == "" || inputFa[0].Name == chrom {
			chromInInput = true
			answer = FaToPfa(inputFa[0], start, end)
		}
	} else {
		if chrom == "" {
			log.Fatalf("Error: expecting a Chrom argument for multifasta input.")
		}

		for _, seq := range inputFa {
			if seq.Name == chrom {
				chromInInput = true
				answer = FaToPfa(seq, start, end)
				break
			}
		}
	}

	if chromInInput == false {
		log.Fatalf("Error: input sequence name does not match requested chrom.")
	}

	return answer
}

// RandSeq returns a randomly-generated pFasta of sequence length 'length' and name 'name' where each base sums to 1 across its probabilities
func RandSeq(length int, name string, seedSet bool, setSeed int64, randSource *rand.Rand) PFasta {
	var source *rand.Rand
	if !seedSet {
		source = rand.New(rand.NewSource(setSeed))
	} else {
		source = randSource
	}
	answer := PFasta{Name: name, Seq: make([]pDna.Float32Base, length)}
	for pos := 0; pos < length; pos++ {
		temp := pDna.RandBase(true, setSeed, source)
		answer.Seq[pos] = temp
	}
	return answer
}
