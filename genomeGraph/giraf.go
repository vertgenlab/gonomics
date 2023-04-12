package genomeGraph

import (
	"fmt"
	"log"
	"math/rand"

	"github.com/vertgenlab/gonomics/cigar"
	"github.com/vertgenlab/gonomics/common"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/giraf"
	"github.com/vertgenlab/gonomics/numbers"
)

func RandGiraf(graph *GenomeGraph, numReads int, readLen int, randSeed int64) []*giraf.Giraf {
	answer := make([]*giraf.Giraf, 0)
	var curr *giraf.Giraf

	var totalBases = BasesInGraph(graph)
	if readLen > totalBases { // not a perfect check, but best we can do without a search algorithm
		log.Fatal("Cannot request more bases than is present in graph")
	}
	rand.Seed(randSeed)

	for len(answer) < numReads {
		nodeIdx, pos := RandLocationFast(graph, totalBases)
		path, endPos, seq := RandPathFwd(graph, nodeIdx, pos, readLen)
		strand := rand.Intn(2) == 0

		if (len(seq) == readLen) && (dna.CountBaseInterval(seq, dna.N, 0, readLen) == 0) {
			girafPath := giraf.Path{
				TStart: int(pos),
				Nodes:  path,
				TEnd:   int(endPos)}

			qual, alnScore, mapQ := generateDiverseQuals(readLen)

			curr = &giraf.Giraf{
				QName:     fmt.Sprintf("%d_%d_%d_%d_%c", path[0], pos+1, path[len(path)-1], endPos+1, common.StrandToRune(strand)),
				QStart:    0,
				QEnd:      readLen,
				PosStrand: strand,
				Path:      girafPath,
				Cigar:     []cigar.ByteCigar{{RunLen: uint16(readLen), Op: 'M'}}, // tmp cigar until giraf cigars have been implemented
				AlnScore:  alnScore,                                              // placeholder
				MapQ:      mapQ,                                                  // placeholder
				Seq:       seq,
				Qual:      qual,
				Notes:     nil}

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
		alnScore = numbers.RandIntInRange(6000, 8000)
		mapQ = uint8(rand.Intn(5))
	case scoreProb < 10: // 10% of bases will have alnScore between 8k-10k, and mapQ below 15
		alnScore = numbers.RandIntInRange(8000, 10000)
		mapQ = uint8(numbers.RandIntInRange(5, 15))
	case scoreProb < 20: // 20% of bases will have alnScore between 10k-15k, and mapQ below 30
		alnScore = numbers.RandIntInRange(10000, 15000)
		mapQ = uint8(numbers.RandIntInRange(15, 30))
	default: // 80% of bases will have alnScore between 15k-20k, and mapQ below 40
		alnScore = numbers.RandIntInRange(15000, 20000)
		mapQ = uint8(numbers.RandIntInRange(30, 40))
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

// TODO: simulate indels
func RandSomaticMutations(graph *GenomeGraph, reads []*giraf.Giraf, numSomaticSNV int, AlleleFrequency float64, randSeed int64) ([]uint32, []uint32) {
	var totalBases = BasesInGraph(graph)
	var mutationNode, mutationPos []uint32
	var nodeIdx uint32
	var pos, readPos uint32
	rand.Seed(randSeed)

	for i := 0; i < numSomaticSNV; i++ {
		nodeIdx, pos = RandLocationFast(graph, totalBases)
		mutationNode = append(mutationNode, nodeIdx)
		mutationPos = append(mutationPos, pos)
		var mutantBase dna.Base = 4
		for j := 0; j < len(reads); j++ {
			for k := 0; k < len(reads[j].Path.Nodes); k++ {
				if reads[j].Path.Nodes[k] == nodeIdx {
					if reads[j].Path.Nodes[0] == nodeIdx && reads[j].Path.TStart > int(pos) {
						continue
					}
					if reads[j].Path.Nodes[len(reads[j].Path.Nodes)-1] == nodeIdx && reads[j].Path.TEnd < int(pos) {
						continue
					}
					readPos = NodePosToReadPos(graph, reads[j], nodeIdx, pos)
					if int(readPos) >= len(reads[j].Seq) {
						continue
					}
					if mutantBase == 4 {
						base := reads[j].Seq[readPos]
						for {
							mutantBase = dna.Base(rand.Intn(4))
							if mutantBase != base {
								break
							}
						}
					}

					randProb := float64(rand.Intn(100)) / 100
					if randProb <= AlleleFrequency {
						reads[j].Seq[readPos] = mutantBase
					}
				}
			}
		}
	}
	return mutationNode, mutationPos
}

func NodePosToReadPos(graph *GenomeGraph, read *giraf.Giraf, node uint32, pos uint32) uint32 {
	var posInPath int
	var readPos uint32 = 0

	for i := 0; i < len(read.Path.Nodes); i++ {
		if read.Path.Nodes[i] == node {
			posInPath = i
			break
		}
	}

	for i := 0; i < posInPath; i++ {
		readPos += uint32(len(graph.Nodes[read.Path.Nodes[i]].Seq))
	}

	readPos += pos
	readPos -= uint32(read.Path.TStart)

	return readPos
}
