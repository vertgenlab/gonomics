package genomeGraph

import (
	"fmt"
	"math/rand"

	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/dna/dnaTwoBit"
	"github.com/vertgenlab/gonomics/fastq"
	"github.com/vertgenlab/gonomics/numbers"
	"github.com/vertgenlab/gonomics/numbers/parse"
)

func RandomPairedReads(genome *GenomeGraph, readLength int, numReads int, numChanges int) []fastq.PairedEnd {
	var answer []fastq.PairedEnd = make([]fastq.PairedEnd, numReads)
	var seq []dna.Base
	var path []uint32
	var nodeIdx, start1, endPos uint32
	var strand bool
	var totalBases = BasesInGraph(genome)
	var fragLen int
	seed := rand.New(rand.NewSource(0))

	for i := 0; i < numReads; {
		strand = numbers.RandIntInRange(0, 2, seed) == 0
		fragLen = numbers.RandIntInRange(300, 500, seed)
		nodeIdx, start1 = RandLocationFast(genome, totalBases, seed)
		path, endPos, seq = RandPathFwd(genome, nodeIdx, start1, fragLen, seed)

		if (len(seq) == fragLen) && (dna.CountBaseInterval(seq, dna.N, 0, readLength) == 0) {
			curr := fastq.PairedEnd{Fwd: fastq.Fastq{}, Rev: fastq.Fastq{}}
			curr.Fwd.Name = fmt.Sprintf("%d_%d_%d_%d_%c_R: 1", path[0], start1, path[len(path)-1], start1+uint32(readLength), parse.StrandToRune(strand))
			curr.Fwd.Seq = make([]dna.Base, readLength)
			copy(curr.Fwd.Seq, seq[:readLength])

			curr.Fwd.Qual = fastq.ToQualUint8(generateDiverseFakeQual(readLength, seed))
			curr.Rev.Name = fmt.Sprintf("%d_%d_%d_%d_%c_R: 2", path[0], start1+uint32(fragLen-readLength), path[len(path)-1], endPos, parse.StrandToRune(strand))
			curr.Rev.Seq = make([]dna.Base, readLength)
			copy(curr.Rev.Seq, seq[uint32(fragLen-readLength):])
			curr.Rev.Qual = fastq.ToQualUint8(generateDiverseFakeQual(readLength, seed))
			if !strand {
				fastq.ReverseComplement(curr.Fwd)
			} else {
				fastq.ReverseComplement(curr.Rev)
			}
			mutate(curr.Fwd.Seq, numChanges, seed)
			mutate(curr.Rev.Seq, numChanges, seed)
			answer[i] = curr
			i++
		}
	}
	return answer
}

func RandLocationFast(genome *GenomeGraph, totalBases int, seed *rand.Rand) (uint32, uint32) {
	var rand int = numbers.RandIntInRange(0, totalBases, seed)
	for i := 0; i < len(genome.Nodes); i++ {
		if rand < genome.Nodes[i].SeqTwoBit.Len {
			return uint32(i), uint32(rand)
		} else {
			rand -= genome.Nodes[i].SeqTwoBit.Len
		}
	}
	return 0, 0 //needed for compiler, should not get here
}

func RandPathFwd(genome *GenomeGraph, nodeIdx uint32, pos uint32, length int, seed *rand.Rand) ([]uint32, uint32, []dna.Base) {
	var answer []dna.Base = make([]dna.Base, 0, length)
	var i int = 0
	for i = 0; i < length && int(pos) < genome.Nodes[nodeIdx].SeqTwoBit.Len; i, pos = i+1, pos+1 {
		answer = append(answer, dnaTwoBit.GetBase(genome.Nodes[nodeIdx].SeqTwoBit, uint(pos)))
	}
	if i == length || len(genome.Nodes[nodeIdx].Next) == 0 {
		return []uint32{nodeIdx}, pos, answer
	} else {
		edgeIdx := numbers.RandIntInRange(0, len(genome.Nodes[nodeIdx].Next), seed)
		return randPathFwdHelper(genome, genome.Nodes[nodeIdx].Next[edgeIdx].Dest.Id, length, answer, []uint32{nodeIdx}, seed)
	}
}

func randPathFwdHelper(genome *GenomeGraph, nodeIdx uint32, length int, progress []dna.Base, path []uint32, seed *rand.Rand) ([]uint32, uint32, []dna.Base) {
	var pos uint32 = 0
	for pos = 0; len(progress) < length && int(pos) < genome.Nodes[nodeIdx].SeqTwoBit.Len; pos = pos + 1 {
		progress = append(progress, dnaTwoBit.GetBase(genome.Nodes[nodeIdx].SeqTwoBit, uint(pos)))
	}
	if len(progress) == length || len(genome.Nodes[nodeIdx].Next) == 0 {
		return append(path, nodeIdx), pos, progress
	} else {
		edgeIdx := numbers.RandIntInRange(0, len(genome.Nodes[nodeIdx].Next), seed)
		return randPathFwdHelper(genome, genome.Nodes[nodeIdx].Next[edgeIdx].Dest.Id, length, progress, append(path, nodeIdx), seed)
	}
}

func RandomReads(genome *GenomeGraph, readLength int, numReads int, numChanges int) []fastq.Fastq {
	var answer []fastq.Fastq = make([]fastq.Fastq, numReads)
	var seq []dna.Base
	var path []uint32
	var nodeIdx, pos, endPos uint32
	var strand bool
	var totalBases = BasesInGraph(genome)
	seed := rand.New(rand.NewSource(0))
	for i := 0; i < numReads; {
		nodeIdx, pos = RandLocationFast(genome, totalBases, seed)
		path, endPos, seq = RandPathFwd(genome, nodeIdx, pos, readLength, seed)
		strand = numbers.RandIntInRange(0, 2, seed) == 0
		if (len(seq) == readLength) && (dna.CountBaseInterval(seq, dna.N, 0, readLength) == 0) {
			curr := fastq.Fastq{}
			curr.Name = fmt.Sprintf("%d_%d_%d_%d_%c", path[0], pos+1, path[len(path)-1], endPos+1, parse.StrandToRune(strand))
			curr.Seq = seq
			curr.Qual = fastq.ToQualUint8(generateDiverseFakeQual(readLength, seed))
			if !strand {
				dna.ReverseComplement(curr.Seq)
			}
			mutate(curr.Seq, numChanges, seed)
			answer[i] = curr
			i++
		}
	}
	return answer
}

func mutate(sequence []dna.Base, numChanges int, seed *rand.Rand) {
	possibleBases := []dna.Base{0, 1, 2, 3}
	for i := 0; i < numChanges; i++ {
		sequence[numbers.RandIntInRange(0, len(sequence), rand.New(rand.NewSource(0)))] = possibleBases[numbers.RandIntInRange(0, len(possibleBases), seed)]
	}
}

func generateDiverseFakeQual(length int, seed *rand.Rand) []rune {
	var answer []rune = make([]rune, length)
	//var asci = []rune{'!', '#', '$', '%', '&', '(', ')', '*', '+', '`', '-', '.', '/', '0', '1', '2', '3', '4', '5', '6', '7', '8', '9', ':', ';', '<', '=', '>', '?', '@', 'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J'}
	var asci = []rune{'F', ',', 'F', ':'}
	for i := 0; i < length; i++ {
		answer[i] = asci[numbers.RandIntInRange(0, len(asci), seed)]
	}
	return answer
}
