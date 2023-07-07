package popgen

import (
	"math"

	"github.com/vertgenlab/gonomics/align"
	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/numbers"
)

// Dunn takes a region of the genome, a multiple alignment, a set of groups, and an option to realign
// the bed region of the alignment before calculating the Dunn Index.
// The returns are the dunn index as a float64, the number of segregating sites considered, and missing group members as a string.
// Mathematical details of the Dunn Index are described at https://en.wikipedia.org/wiki/Dunn_index.
func Dunn(b bed.Bed, aln []fasta.Fasta, g []*Group, realign bool) (float64, int, string) {
	var maxIntra int = 0
	var minInter int = 0
	var missing string = ""
	var tmp2Fa, tmp3Fa []fasta.Fasta
	//bLen := b.ChromEnd - b.ChromStart

	//generate subFa - filters out alignment positions with gaps or complete identity
	alnPos := fasta.RefPosToAlnPos(aln[0], b.ChromStart)
	alnEnd := fasta.RefPosToAlnPos(aln[0], b.ChromEnd) // could be gaps between ChromStart and ChromEnd
	tmpFa := fasta.CopySubset(aln, alnPos, alnEnd)
	if realign {
		tmp2Fa = fasta.RemoveGaps(tmpFa)
		//fasta.AllToUpper(tmp2Fa)
		tmp2Fa = FilterMultByGroup(tmp2Fa, g)
		tmp3Fa = align.AllSeqAffine(tmp2Fa, align.DefaultScoreMatrix, -400, -30)
		//tmp3Fa = FilterMultByGroup(tmp3Fa, g)
	} else {
		tmp2Fa = fasta.RemoveMissingMult(tmpFa)
		tmp3Fa = FilterMultByGroup(tmp2Fa, g)
	}
	if len(tmp3Fa) == 0 {
		return -1.0, 0, missing //dunn could not be calculated
	}
	subFa := fasta.DistColumn(tmp3Fa)

	missing = FindMissingGroupMembers(subFa, g)

	for i := 0; i < len(g); i++ {
		maxIntra = numbers.Max(maxIntra, findMaxIntra(subFa, g[i]))
	}

	minInter = findMinInter(g, subFa)
	//fmt.Printf("%d\t%d\t%d\t%d\t%g\n", b.ChromStart, fasta.NumSegregatingSites(subFa), minInter, maxIntra, float64(minInter)/float64(maxIntra))
	return float64(minInter) / float64(maxIntra), fasta.NumSegregatingSites(subFa), missing
}

// findMaxIntra is a helper function of Dunn that calculates the Max pairwise sequence distance between two sequences of a multiFa alignment that are part of the same Group.
func findMaxIntra(subFa []fasta.Fasta, g *Group) int {
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

// findMinInter is a helper function of Dunn that calculates the minimum pairwise sequence distance between two sequences of a multiFa alignment that are part of different groups.
func findMinInter(g []*Group, subFa []fasta.Fasta) int {
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
