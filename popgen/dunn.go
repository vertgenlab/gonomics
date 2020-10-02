package popgen

import (
	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/numbers"
	//"log"
	"math"
)

//Dunn returns the dunn index as a float64 and missing group members as a string.
//Mathematical details of the Dunn Index are described at https://en.wikipedia.org/wiki/Dunn_index.
func Dunn(b *bed.Bed, aln []*fasta.Fasta, g []*Group) (float64, string) {
	var maxIntra int = 0
	var minInter int = 0
	var missing string = ""
	bLen := b.ChromEnd - b.ChromStart

	//generate subFa - filters out alignment positions with gaps or complete identity
	alnPos := fasta.RefPosToAlnPos(aln[0], int(b.ChromStart))
	//fmt.Printf("RePos Done.\n")
	tmpFa := fasta.CopySubset(aln, alnPos, alnPos+int(bLen))
	tmp2Fa := fasta.RemoveMissingMult(tmpFa)
	tmp3Fa := FilterMultByGroup(tmp2Fa, g)
	//fmt.Printf("Filter range ok.\n")
	subFa := fasta.DistColumn(tmp3Fa)
	//fmt.Printf("DistBase completed.\n")

	missing = FindMissingGroupMembers(subFa, g)

	for i := 0; i < len(g); i++ {
		maxIntra = numbers.Max(maxIntra, FindMaxIntra(subFa, g[i], b))
	}

	minInter = FindMinInter(g, subFa)

	return (float64(minInter) / float64(maxIntra)), missing
}

//FindMaxIntra is a helper function of Dunn that calculates the Max pairwise sequence distance between two sequences of a multiFa alignment that are part of the same Group.
func FindMaxIntra(subFa []*fasta.Fasta, g *Group, b *bed.Bed) int {
	var answer int = 0
	var group1index int
	var group2index int
	for i := 0; i < len(g.Members); i++ {
		for j := i + 1; j < len(g.Members); j++ {
			group1index = fasta.FindFaIndex(subFa, g.Members[i])
			group2index = fasta.FindFaIndex(subFa, g.Members[j])
			if group1index != -1 && group2index != -1 {
				answer = numbers.Max(answer, dna.Dist(subFa[group1index].Seq, subFa[group2index].Seq))
			}
		}
	}
	return answer
}

//FindMinInter is a helper function of Dunn that calculates the minimum pairwise sequence distance between two sequences of a multiFa alignment that are part of different groups.
func FindMinInter(g []*Group, subFa []*fasta.Fasta) int {
	var answer int = math.MaxInt64
	var group1index, group2index int
	for i := 0; i < len(g[0].Members); i++ {
		for j := 0; j < len(g[1].Members); j++ {
			group1index = fasta.FindFaIndex(subFa, g[0].Members[i])
			group2index = fasta.FindFaIndex(subFa, g[1].Members[j])
			if group1index != -1 && group2index != -1 {
				answer = numbers.Min(answer, dna.Dist(subFa[group1index].Seq, subFa[group2index].Seq))
			}
		}
	}
	return int(answer)
}
