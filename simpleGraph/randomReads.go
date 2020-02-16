package simpleGraph

import (
	"fmt"
	"github.com/vertgenlab/gonomics/common"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/fastq"
	"log"
	"math/rand"
)

func PairedEndRandomReads(genome []*Node, readLength int, numReads int, numChanges int) []*fastq.PairedEnd {
	simReads := RandomReads(genome, readLength, numReads, numChanges)
	log.Printf("length of single generated reads: %d\n", len(simReads))
	simPairs := make([]*fastq.PairedEnd, len(simReads))
	for i := 0; i < len(simReads); i++ {
		curr := fastq.PairedEnd{}
		curr = singleToPair(simReads[i])
		simPairs[i] = &curr
	}
	return simPairs
}

func singleToPair(inFq *fastq.Fastq) fastq.PairedEnd {
	fwd := &fastq.Fastq{Name: fmt.Sprintf("%s 1", inFq.Name), Seq: make([]dna.Base, len(inFq.Seq)), Qual: make([]rune, len(inFq.Qual))}
	copy(fwd.Seq, inFq.Seq)
	copy(fwd.Qual, inFq.Qual)
	rev := &fastq.Fastq{Name: fmt.Sprintf("%s 2", inFq.Name), Seq: make([]dna.Base, len(inFq.Seq)), Qual: make([]rune, len(inFq.Qual))}
	copy(rev.Seq, inFq.Seq)
	copy(rev.Qual, inFq.Qual)
	fastq.ReverseComplement(rev)

	fqPair := fastq.PairedEnd{Fwd: nil, Rev: nil}
	fqPair.Fwd = fwd
	fqPair.Rev = rev
	return fqPair
}

func mutateSingleRead(read fastq.Fastq, location int64, size int) fastq.Fastq {
	var editGenome int = randIntInRange(0, 3)
	possibleBases := []dna.Base{dna.A, dna.C, dna.G, dna.T}
	extra := make([]dna.Base, size)
	if editGenome == 0 {
		mutatePos(read.Seq, int(location))
	} else if editGenome == 1 {
		for base := 0; base < size; base++ {
			extra[base] = possibleBases[randIntInRange(0, len(possibleBases))]
		}
		dna.Insert(read.Seq, location, extra)

	} else {
		dna.Delete(read.Seq, location, location+int64(size))
	}
	return read
}

func MutateRandomReads(genome []*Node, readLength int, numReads int, numChanges int) []*fastq.Fastq {
	var answer []*fastq.Fastq = make([]*fastq.Fastq, numReads)
	var start int
	var chromIdx int
	var strand bool
	var i int

	for i = 0; i < numReads; {
		chromIdx = randIntInRange(0, len(genome))
		size := randIntInRange(2, 8)
		location := int64(randIntInRange(0, readLength-size))
		start = randIntInRange(0, len(genome[chromIdx].Seq)-readLength)
		strand = randIntInRange(0, 2) == 0
		if dna.CountBaseInterval(genome[chromIdx].Seq, dna.N, start, start+readLength) == 0 {
			curr := fastq.Fastq{}
			curr.Name = fmt.Sprintf("%d_%d_%d_%c", genome[chromIdx].Id, start, start+readLength, common.StrandToRune(strand))
			curr.Seq = make([]dna.Base, readLength)
			copy(curr.Seq, genome[chromIdx].Seq[start:start+readLength])
			curr.Qual = generateDiverseFakeQual(readLength)
			curr = mutateSingleRead(curr, location, size)
			if !strand {
				dna.ReverseComplement(curr.Seq)
			}

			answer[i] = &curr
			i++
		}
	}
	return answer
}

func RandomReads(genome []*Node, readLength int, numReads int, numChanges int) []*fastq.Fastq {
	var answer []*fastq.Fastq = make([]*fastq.Fastq, numReads)
	var start int
	var chromIdx int
	var strand bool
	for i := 0; i < numReads; {
		chromIdx = randIntInRange(0, len(genome))
		start = randIntInRange(0, len(genome[chromIdx].Seq)-readLength)
		strand = randIntInRange(0, 2) == 0
		if dna.CountBaseInterval(genome[chromIdx].Seq, dna.N, start, start+readLength) == 0 {
			curr := fastq.Fastq{}
			curr.Name = fmt.Sprintf("%d_%d_%d_%c", genome[chromIdx].Id, start, start+readLength, common.StrandToRune(strand))
			curr.Seq = make([]dna.Base, readLength)
			copy(curr.Seq, genome[chromIdx].Seq[start:start+readLength])
			curr.Qual = generateDiverseFakeQual(readLength)
			if !strand {
				dna.ReverseComplement(curr.Seq)
			}
			mutate(curr.Seq, numChanges)
			answer[i] = &curr
			i++
		}
	}
	return answer
}

func RandomFastqGen(genome []*fasta.Fasta, readLength int, numReads int) []*fastq.Fastq {
	var answer []*fastq.Fastq = make([]*fastq.Fastq, numReads)
	var start int
	var chromIdx int
	var readName string
	var strand bool
	var qual []rune
	var seq []dna.Base
	for i := 0; i < numReads; {
		chromIdx = randIntInRange(0, len(genome))
		start = randIntInRange(0, len(genome[chromIdx].Seq)-readLength)
		strand = randIntInRange(0, 2) == 0
		if dna.CountBaseInterval(genome[chromIdx].Seq, dna.N, start, start+readLength) == 0 {
			readName = fmt.Sprintf("%s_%d_%d_%c", genome[chromIdx].Name, start, start+readLength, common.StrandToRune(strand))
			qual = generateDiverseFakeQual(readLength)
			seq = genome[chromIdx].Seq[start : start+readLength]
			if !strand {
				dna.ReverseComplement(seq)
			}
			answer[i] = &fastq.Fastq{Name: readName, Seq: seq, Qual: qual}
			i++
		}
	}
	return answer
}

func randIntInRange(x int, y int) int {
	return int(rand.Float64()*float64(y-x)) + x
}

func mutate(sequence []dna.Base, numChanges int) {
	possibleBases := []dna.Base{0, 1, 2, 3}
	for i := 0; i < numChanges; i++ {
		sequence[randIntInRange(0, len(sequence))] = possibleBases[randIntInRange(0, len(possibleBases))]
	}
}

func mutatePos(seq []dna.Base, pos int) {
	possibleBases := []dna.Base{dna.A, dna.C, dna.G, dna.T}
	newBase := possibleBases[randIntInRange(0, len(possibleBases))]
	if newBase == seq[pos] {
		mutatePos(seq, pos)
	} else {
		seq[pos] = newBase
	}
}

func RandomReadOneMutation(genome []*Node, readLength int, mutantPos int) *fastq.Fastq {
	var answer *fastq.Fastq = nil
	var start int
	var chromIdx int
	var strand bool
	var readOk bool = false

	for !readOk {
		chromIdx = randIntInRange(0, len(genome))
		start = randIntInRange(0, len(genome[chromIdx].Seq)-readLength)
		if dna.CountBaseInterval(genome[chromIdx].Seq, dna.N, start, start+readLength) == 0 {
			strand = randIntInRange(0, 2) == 0
			curr := fastq.Fastq{}
			curr.Name = fmt.Sprintf("%d_%d_%d_%c", genome[chromIdx].Id, start, start+readLength, common.StrandToRune(strand))
			curr.Seq = make([]dna.Base, readLength)
			copy(curr.Seq, genome[chromIdx].Seq[start:start+readLength])
			curr.Qual = generateFakeQual(readLength)
			if !strand {
				dna.ReverseComplement(curr.Seq)
			}
			mutatePos(curr.Seq, mutantPos)
			answer = &curr
			readOk = true
		}
	}
	return answer
}

func generateDiverseFakeQual(length int) []rune {
	var answer []rune = make([]rune, length)
	//var asci = []rune{'!', '#', '$', '%', '&', '(', ')', '*', '+', '`', '-', '.', '/', '0', '1', '2', '3', '4', '5', '6', '7', '8', '9', ':', ';', '<', '=', '>', '?', '@', 'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J'}
	var asci = []rune{'F', ',', 'F', ':'}
	for i := 0; i < length; i++ {
		answer[i] = asci[randIntInRange(0, len(asci))]
	}
	return answer
}

func generateFakeQual(length int) []rune {
	var answer []rune = make([]rune, length)
	for i := 0; i < length; i++ {
		answer[i] = 'J'
	}
	return answer
}

func generateQual(bases []dna.Base) []rune {
	var ans []rune
	for i := 0; i < len(bases); i++ {
		ans = append(ans, 'J')
	}
	return ans
}
