package simulate

import (
	"fmt"
	"github.com/vertgenlab/gonomics/cigar"
	"github.com/vertgenlab/gonomics/common"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/giraf"
	"github.com/vertgenlab/gonomics/numbers"
	"github.com/vertgenlab/gonomics/simpleGraph"
	"log"
	"math/rand"
)

func RandGiraf(graph *simpleGraph.SimpleGraph, numReads int, readLen int, randSeed int64) []*giraf.Giraf {
	answer := make([]*giraf.Giraf, 0)
	var curr *giraf.Giraf

	var totalBases = ByteBasesInGraph(graph)
	rand.Seed(randSeed)


	for ; len(answer) < numReads ; {
		nodeIdx, pos := ByteRandLocationFast(graph, totalBases)
		path, endPos, seq := ByteRandPathFwd(graph, nodeIdx, pos, readLen)
		strand := rand.Intn(2) == 0

		if (len(seq) == readLen) && (dna.CountBaseInterval(seq, dna.N, 0, readLen) == 0) {
			girafPath := &giraf.Path{
				TStart: int(pos+1),
				Nodes: path,
				TEnd: int(endPos)}

			qual, alnScore, mapQ := generateDiverseQuals(readLen)

			curr = &giraf.Giraf{
				QName: fmt.Sprintf("%d_%d_%d_%d_%c", path[0], pos+1, path[len(path)-1], endPos+1, common.StrandToRune(strand)),
				QStart: 0,
				QEnd: readLen,
				PosStrand: strand,
				Path: girafPath,
				Aln: []*cigar.Cigar{{int64(readLen), 'M'}}, // tmp cigar until giraf cigars have been implemented
				AlnScore: alnScore, // placeholder
				MapQ: mapQ, // placeholder
				Seq: seq,
				Qual: qual,
				Notes: nil}


			if !strand {
				dna.ReverseComplement(curr.Seq)
			}
			answer = append(answer, curr)
		}
	}
	return answer
}

func generateDiverseQuals(readLen int) ([]uint8, int, uint8) {
	answer := make([]uint8, readLen)
	var alnScore int
	var mapQ uint8

	scoreProb := rand.Intn(100)
	switch {
	case scoreProb == 0: // 1% of bases will have alnScore between 6k-8k, and mapQ below 5
		alnScore = 6000 + rand.Intn(2000)
		mapQ = uint8(rand.Intn(5))
	case scoreProb < 10: // 10% of bases will have alnScore between 8k-10k, and mapQ below 15
		alnScore = 8000 + rand.Intn(2000)
		mapQ = 5 + uint8(rand.Intn(10))
	case scoreProb < 20: // 20% of bases will have alnScore between 10k-15k, and mapQ below 30
		alnScore = 10000 + rand.Intn(5000)
		mapQ = 15 + uint8(rand.Intn(15))
	default: // 80% of bases will have alnScore between 15k-20k, and mapQ below 40
		alnScore = 15000 + rand.Intn(5000)
		mapQ = 30 + uint8(rand.Intn(10))
	}

	for i := 0; i < readLen; i++ {
		scoreProb := rand.Intn(100)
		scoreBase := rand.Intn(10)

		switch {
		case scoreProb == 0: // 1% of bases will have qual < 10
			answer[i] = uint8(scoreBase)
		case scoreProb < 10: // 10% of bases will have qual < 20
			answer[i] = uint8(scoreBase + 10)
		case scoreProb < 20: // 20% of bases will have qual < 30
			answer[i] = uint8(scoreBase + 20)
		default: // 80% of bases will have qual > 30
			answer[i] = uint8(scoreBase + 30)
		}
	}
	return answer, alnScore, mapQ
}

func RandSomaticMutations(graph *simpleGraph.SimpleGraph, reads []*giraf.Giraf, numSomaticSNV int, numSomaticIns int, numSomaticDel int, randSeed int64) []*giraf.Giraf {
	//var totalBases = ByteBasesInGraph(graph)
	rand.Seed(randSeed)

	return reads
}

func ByteBasesInGraph(g *simpleGraph.SimpleGraph) int {
	var i, baseCount int = 0, 0
	for i = 0; i < len(g.Nodes); i++ {
		baseCount += len(g.Nodes[i].Seq)
	}
	return baseCount
}

func ByteRandLocationFast(genome *simpleGraph.SimpleGraph, totalBases int) (uint32, uint32) {
	var rand int = numbers.RandIntInRange(0, totalBases)
	for i := 0; i < len(genome.Nodes); i++ {
		if rand < len(genome.Nodes[i].Seq) {
			return uint32(i), uint32(rand)
		} else {
			rand -= len(genome.Nodes[i].Seq)
		}
	}
	log.Fatal("Error: trouble selecting a random location in the graph\n")
	return 0, 0 //needed for compiler, should not get here
}

func ByteRandPathFwd(genome *simpleGraph.SimpleGraph, nodeIdx uint32, pos uint32, length int) ([]uint32, uint32, []dna.Base) {
	var answer []dna.Base = make([]dna.Base, 0, length)
	var i int = 0
	for i = 0; i < length && int(pos) < len(genome.Nodes[nodeIdx].Seq); i, pos = i+1, pos+1 {
		answer = append(answer, genome.Nodes[nodeIdx].Seq[pos])
	}
	if i == length || len(genome.Nodes[nodeIdx].Next) == 0 {
		return []uint32{nodeIdx}, pos, answer
	} else {
		edgeIdx := numbers.RandIntInRange(0, len(genome.Nodes[nodeIdx].Next))
		return randPathFwdHelper(genome, genome.Nodes[nodeIdx].Next[edgeIdx].Dest.Id, length, answer, []uint32{nodeIdx})
	}
}

func randPathFwdHelper(genome *simpleGraph.SimpleGraph, nodeIdx uint32, length int, progress []dna.Base, path []uint32) ([]uint32, uint32, []dna.Base) {
	var pos uint32 = 0
	for pos = 0; len(progress) < length && int(pos) < len(genome.Nodes[nodeIdx].Seq); pos = pos + 1 {
		progress = append(progress, genome.Nodes[nodeIdx].Seq[pos])
	}
	if len(progress) == length || len(genome.Nodes[nodeIdx].Next) == 0 {
		return append(path, nodeIdx), pos, progress
	} else {
		edgeIdx := numbers.RandIntInRange(0, len(genome.Nodes[nodeIdx].Next))
		return randPathFwdHelper(genome, genome.Nodes[nodeIdx].Next[edgeIdx].Dest.Id, length, progress, append(path, nodeIdx))
	}
}