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
func Dunn(b *bed.Bed, aln []fasta.Fasta, g []*Group) (float64, string) {
	var maxIntra int = 0
	var minInter int = 0
	var missing string = ""
	bLen := b.ChromEnd - b.ChromStart

	//generate subFa - filters out alignment positions with gaps or complete identity
	alnPos := fasta.RefPosToAlnPos(aln[0], b.ChromStart)
	tmpFa := fasta.CopySubset(aln, alnPos, alnPos+bLen)
	tmp2Fa := fasta.RemoveMissingMult(tmpFa)
	tmp3Fa := FilterMultByGroup(tmp2Fa, g)
	if len(tmp3Fa) == 0 {
		return -1.0, missing//dunn could not be calculated
	}
	subFa := fasta.DistColumn(tmp3Fa)

	missing = FindMissingGroupMembers(subFa, g)

	for i := 0; i < len(g); i++ {
		maxIntra = numbers.Max(maxIntra, FindMaxIntra(subFa, g[i]))
	}

	minInter = FindMinInter(g, subFa)

	return float64(minInter) / float64(maxIntra), missing
}

//FindMaxIntra is a helper function of Dunn that calculates the Max pairwise sequence distance between two sequences of a multiFa alignment that are part of the same Group.
func FindMaxIntra(subFa []fasta.Fasta, g *Group) int {
	var answer int = 0
	var faI, faJ []dna.Base
	var faFoundI, faFoundJ bool
	faMap := fasta.ToMap(subFa)
	for i := 0; i < len(g.Members); i++ {
		for j := i + 1; j < len(g.Members); j++ {
			faI, faFoundI = faMap[g.Members[i]]
			faJ, faFoundJ = faMap[g.Members[j]]
			if faFoundI && faFoundJ {
				answer = numbers.Max(answer, dna.Dist(faI, faJ))
			}
		}
	}
	return answer
}

//FindMinInter is a helper function of Dunn that calculates the minimum pairwise sequence distance between two sequences of a multiFa alignment that are part of different groups.
func FindMinInter(g []*Group, subFa []fasta.Fasta) int {
	var answer int = math.MaxInt64
	faMap := fasta.ToMap(subFa)
	var faI, faJ []dna.Base
	var faFoundI, faFoundJ bool
	for i := 0; i < len(g[0].Members); i++ {
		for j := 0; j < len(g[1].Members); j++ {
			faI, faFoundI = faMap[g[0].Members[i]]
			faJ, faFoundJ = faMap[g[1].Members[j]]
			if faFoundI && faFoundJ {
				answer = numbers.Min(answer, dna.Dist(faI, faJ))
			}
		}
	}
	return answer
}
