package simpleGraph

import (
	"fmt"
	"github.com/vertgenlab/gonomics/common"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/dnaTwoBit"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/fastq"
	"github.com/vertgenlab/gonomics/numbers"
	"log"
	"math/rand"
)

//TODO: Check for bugs, failed align tests so it's unclear if this function is wrong or alignment function is wrong...
func GraphRandomReads(genome []*Node, readLength int, numReads int, numChanges int) []*fastq.Fastq {
	var answer []*fastq.Fastq = make([]*fastq.Fastq, numReads)
	var start int
	var chromIdx int
	var strand bool
	for i := 0; i < numReads; {
		strand = randIntInRange(0, 2) == 0
		chromIdx = randIntInRange(0, len(genome))
		//curr.Name = fmt.Sprintf("%d_%d_%d_%c", genome[chromIdx].Id, start+1, start+1+readLength, common.StrandToRune(strand))
		if dna.CountBaseInterval(genome[chromIdx].Seq, dna.N, start, common.Min(start+readLength, len(genome[chromIdx].Seq))) == 0 {
			curr := fastq.Fastq{}
			curr.Seq = make([]dna.Base, readLength)
			if len(genome[chromIdx].Seq) < readLength {
				start = randIntInRange(0, len(genome[chromIdx].Seq))
				curr.Name = fmt.Sprintf("%d_%d_%d_%c", genome[chromIdx].Id, start+1, len(genome[chromIdx].Seq), common.StrandToRune(strand))
				curr = getSeqRandHelp(genome[chromIdx], []dna.Base{}, start, readLength, curr)
			} else {
				start = randIntInRange(0, len(genome[chromIdx].Seq)-readLength)
				curr.Name = fmt.Sprintf("%d_%d_%d_%c", genome[chromIdx].Id, start+1, start+1+readLength, common.StrandToRune(strand))
				copy(curr.Seq, genome[chromIdx].Seq[start:start+readLength])
			}
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

func getSeqRandHelp(curr *Node, seq []dna.Base, start int, readLength int, answer fastq.Fastq) fastq.Fastq {
	//var answer *fastq.Fastq
	var availableBases int = len(curr.Seq) - start + len(seq)
	var targetLength int = common.Min(availableBases, readLength)
	var basesToTake int = targetLength - len(seq)
	//log.Printf("len(seq)=%d, len(n.Seq)=%d, start=%d, targetLength=%d, basesToTake=%d\n", len(seq), len(curr.Seq), start, targetLength, basesToTake)
	var s []dna.Base = make([]dna.Base, targetLength)
	copy(s[0:len(seq)], seq)
	copy(s[len(seq):targetLength], curr.Seq[start:start+basesToTake])
	if availableBases >= readLength || len(curr.Next) == 0 {
		if dna.CountBaseInterval(s, dna.N, 0, len(s)) == 0 {
			copy(answer.Seq, s)
		}
		return answer
	} else {
		edgeIdx := randIntInRange(0, len(curr.Next))
		answer.Name += fmt.Sprintf(":%d", curr.Next[edgeIdx].Dest.Id)
		answer = getSeqRandHelp(curr.Next[edgeIdx].Dest, s, 0, readLength, answer)
		return answer
	}
}

func RandomPairedReads(genome []*Node, readLength int, numReads int, numChanges int) []*fastq.PairedEnd {
	var answer []*fastq.PairedEnd = make([]*fastq.PairedEnd, numReads)
	var start, start2 int
	var fragLen int
	var chromIdx int
	var strand bool
	for i := 0; i < numReads; {
		chromIdx = randIntInRange(0, len(genome))
		fragLen = randIntInRange(300, 500)
		start = randIntInRange(0, len(genome[chromIdx].Seq)-fragLen)
		start2 = start + fragLen - readLength
		strand = randIntInRange(0, 2) == 0
		if dna.CountBaseInterval(genome[chromIdx].Seq, dna.N, start, start+readLength) == 0 {
			curr := fastq.PairedEnd{Fwd: &fastq.Fastq{}, Rev: &fastq.Fastq{}}

			curr.Fwd.Name = fmt.Sprintf("%d_%d_%d_%c_%d_R: 1", genome[chromIdx].Id, start+1, start+1+readLength, common.StrandToRune(strand), fragLen)
			curr.Fwd.Seq = make([]dna.Base, readLength)
			copy(curr.Fwd.Seq, genome[chromIdx].Seq[start:start+readLength])
			curr.Fwd.Qual = generateDiverseFakeQual(readLength)

			curr.Rev.Name = fmt.Sprintf("%d_%d_%d_%c_%d_R: 2", genome[chromIdx].Id, start2+1, fragLen+start+1, common.StrandToRune(strand), fragLen)
			curr.Rev.Seq = make([]dna.Base, readLength)
			copy(curr.Rev.Seq, genome[chromIdx].Seq[start2:start2+readLength])
			curr.Rev.Qual = generateDiverseFakeQual(readLength)

			if !strand {
				fastq.ReverseComplement(curr.Fwd)
			} else {
				fastq.ReverseComplement(curr.Rev)
			}

			mutate(curr.Fwd.Seq, numChanges)
			mutate(curr.Rev.Seq, numChanges)
			answer[i] = &curr
			i++
		}
	}
	return answer
}

func PairedEndRandomReads(genome *SimpleGraph, readLength int, numReads int, numChanges int) []*fastq.PairedEnd {
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

func RandLocation(genome *SimpleGraph) (uint32, uint32) {
	var totalBases int = BasesInGraph(genome)
	return RandLocationFast(genome, totalBases)
}

func RandLocationFast(genome *SimpleGraph, totalBases int) (uint32, uint32) {
	var rand int = numbers.RandIntInRange(0, totalBases)
	for i := 0; i < len(genome.Nodes); i++ {
		if rand < genome.Nodes[i].SeqTwoBit.Len {
			return uint32(i), uint32(rand)
		} else {
			rand -= genome.Nodes[i].SeqTwoBit.Len
		}
	}
	log.Fatal("Error: trouble selecting a random location in the graph\n")
	return 0, 0 //needed for compiler, should not get here
}

func RandPathFwd(genome *SimpleGraph, nodeIdx uint32, pos uint32, length int) ([]uint32, uint32, []dna.Base) {
	var answer []dna.Base = make([]dna.Base, 0, length)
	var i int = 0
	for i = 0; i < length && int(pos) < genome.Nodes[nodeIdx].SeqTwoBit.Len; i, pos = i+1, pos+1 {
		answer = append(answer, dnaTwoBit.GetBase(genome.Nodes[nodeIdx].SeqTwoBit, uint(pos)))
	}
	if i == length || len(genome.Nodes[nodeIdx].Next) == 0 {
		return []uint32{nodeIdx}, pos, answer
	} else {
		edgeIdx := numbers.RandIntInRange(0, len(genome.Nodes[nodeIdx].Next))
		return randPathFwdHelper(genome, genome.Nodes[nodeIdx].Next[edgeIdx].Dest.Id, length, answer, []uint32{nodeIdx})
	}
}

func randPathFwdHelper(genome *SimpleGraph, nodeIdx uint32, length int, progress []dna.Base, path []uint32) ([]uint32, uint32, []dna.Base) {
	var pos uint32 = 0
	for pos = 0; len(progress) < length && int(pos) < genome.Nodes[nodeIdx].SeqTwoBit.Len; pos = pos + 1 {
		progress = append(progress, dnaTwoBit.GetBase(genome.Nodes[nodeIdx].SeqTwoBit, uint(pos)))
	}
	if len(progress) == length || len(genome.Nodes[nodeIdx].Next) == 0 {
		return append(path, nodeIdx), pos, progress
	} else {
		edgeIdx := numbers.RandIntInRange(0, len(genome.Nodes[nodeIdx].Next))
		return randPathFwdHelper(genome, genome.Nodes[nodeIdx].Next[edgeIdx].Dest.Id, length, progress, append(path, nodeIdx))
	}
}

func RandomReads(genome *SimpleGraph, readLength int, numReads int, numChanges int) []*fastq.Fastq {
	var answer []*fastq.Fastq = make([]*fastq.Fastq, numReads)
	var seq []dna.Base
	var path []uint32
	var nodeIdx, pos, endPos uint32
	var strand bool
	var totalBases = BasesInGraph(genome)
	for i := 0; i < numReads; {
		nodeIdx, pos = RandLocationFast(genome, totalBases)
		path, endPos, seq = RandPathFwd(genome, nodeIdx, pos, readLength)
		strand = numbers.RandIntInRange(0, 2) == 0
		if (len(seq) == readLength) && (dna.CountBaseInterval(seq, dna.N, 0, readLength) == 0) {
			curr := fastq.Fastq{}
			curr.Name = fmt.Sprintf("%d_%d_%d_%d_%c", path[0], pos+1, path[len(path)-1], endPos+1, common.StrandToRune(strand))
			curr.Seq = seq
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
