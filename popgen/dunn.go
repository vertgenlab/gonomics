package popgen

import (
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/common"
	"math"
	"log"
	"github.com/vertgenlab/gonomics/dna"
)

func Dunn(b *bed.Bed, aln []*fasta.Fasta, g []*Group) float64 {
	var maxIntra int = 0
	var minInter int = 0
	bLen := b.ChromEnd - b.ChromStart 

	//generate subFa - filters out alignment positions with gaps or complete identity
	alnPos := fasta.RefPosToAlnPos(aln[0], int(b.ChromStart))
	//fmt.Printf("RePos Done.\n")
	tmpFa := fasta.CopySubset(aln, alnPos, alnPos + int(bLen))
	tmp2Fa := fasta.RemoveMissingMult(tmpFa)
	//fmt.Printf("Filter range ok.\n")
	subFa := fasta.DistColumn(tmp2Fa)
	//fmt.Printf("DistBase completed.\n")

	for i := 0; i < len(g); i++ {
		maxIntra = common.Max(maxIntra, FindMaxIntra(subFa, g[i]))
	}

	minInter = FindMinInter(g, subFa)

	return (float64(minInter) / float64(maxIntra))
}

func FindMaxIntra(subFa []*fasta.Fasta, g *Group) int {
	var answer int = 0
	var group1index int
	var group2index int
	for i := 0; i < len(g.Members); i++ {
		for j := i + 1; j < len(g.Members); j++ { 
			group1index = findFaIndex(subFa, g.Members[i])
			group2index = findFaIndex(subFa, g.Members[j])
			answer = common.Max(answer, dna.Dist(subFa[group1index].Seq, subFa[group2index].Seq))
		}
	}
	return answer
}

func findFaIndex(subFa []*fasta.Fasta, n string) int {
	for i := 0; i < len(subFa); i++ {
		if subFa[i].Name == n {
			return i
		}
	}
	log.Fatalf("Group member: %s not found in alignment.\n", n)
	return -1
}

func FindMinInter(g []*Group, subFa []*fasta.Fasta) int {
	var answer int = math.MaxInt64
	var spec1index, spec2index int
	for i := 0; i < len(g[0].Members); i++ {
		for j := 0; j < len(g[1].Members); j++ {
			spec1index = findFaIndex(subFa, g[0].Members[i])
			spec2index = findFaIndex(subFa, g[1].Members[j])
			answer = common.Min(answer, dna.Dist(subFa[spec1index].Seq, subFa[spec2index].Seq))
		}
	}
	return int(answer)
}
