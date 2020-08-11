package popgen

import (
	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/fasta"
	"math"
)

//TajimaFromBed calculates Tajima's D from a multiFa alignment at a position set specific by an input bed entry.
//Requires a gapless alignment.
func TajimaFromBedNoGroup(b *bed.Bed, aln []*fasta.Fasta) float64 {
	bLen := b.ChromEnd - b.ChromStart
	alnPos := fasta.RefPosToAlnPos(aln[0], int(b.ChromStart))
	tmpFa := fasta.CopySubset(aln, alnPos, alnPos+int(bLen))
	tmpFa = fasta.RemoveMissingMult(tmpFa)

	return Tajima(tmpFa)
}

//TajimaFromBed caculates Tajima's D from a multiFa alignment at a position set specified by an input bed entry.
//This version considers only the alignment entries that are represented in an input Group slice. A second return from this function lists members of the input Group slice not found in the multiFa alignment.
//Requires a gapless alignment.
func TajimaFromBed(b *bed.Bed, aln []*fasta.Fasta, g []*Group) (float64, string) {
	bLen := b.ChromEnd - b.ChromStart
	alnPos := fasta.RefPosToAlnPos(aln[0], int(b.ChromStart))
	tmpFa := fasta.CopySubset(aln, alnPos, alnPos+int(bLen))
	tmp2Fa := fasta.RemoveMissingMult(tmpFa)

	tmp3Fa := FilterMultByGroup(tmp2Fa, g)
	missing := FindMissingGroupMembers(tmp3Fa, g)

	return Tajima(tmpFa), missing
}

//Tajima calculates Tajima's D from an input multiFa alignment block.
func Tajima(aln []*fasta.Fasta) float64 {
	k := calculateTajimaK(aln)
	SoverA := calculateSoverA(aln)
	D := (k - SoverA) / calculateDenominator(aln)
	return D
}

//calculateDenominator is a helper function of Tajima that calculates the denominator of the expression.
func calculateDenominator(aln []*fasta.Fasta) float64 {
	n := float64(len(aln))
	b1 := (n + 1) / (3 * (n - 1))
	b2 := 2 * (math.Pow(n, 2) + n + 3) / (9 * n * (n - 1))
	a1 := calculateA1(aln)
	a2 := calculateA2(aln)

	c1 := b1 - (1 / a1)
	e1 := c1 / a1
	c2 := b2 - (n+2)/(a1*n) + a2/math.Pow(a1, 2)
	e2 := c2 / (math.Pow(a1, 2) + a2)

	S := float64(calculateS(aln))
	return math.Sqrt(e1*S + e2*S*(S-1))
}

//calculateA2 is a helper function of calculateDenominator that returns the value of A2, one parameter of Tajima's D.
func calculateA2(aln []*fasta.Fasta) float64 {
	var a2 float64

	for i := 1; i < (len(aln) - 1); i++ {
		a2 += (1.0 / math.Pow(float64(i), 2))
	}

	return a2
}

//calculateA1 is a helper function of calculateDenominator that returns the value of the A1 parameter.
func calculateA1(aln []*fasta.Fasta) float64 {
	var a1 float64

	for i := 1; i < len(aln)-1; i++ {
		a1 += (1.0 / float64(i))
	}
	return a1
}

//calculateSoverA is a helper function of Tajima that returns the SoverA approximation of theta.
func calculateSoverA(aln []*fasta.Fasta) float64 {
	S := calculateS(aln)
	a1 := calculateA1(aln)
	return (float64(S) / a1)
}

//calculateS returns the number of segregating sites from a multiFa alignment block.
func calculateS(aln []*fasta.Fasta) int {
	S := 0
	var diff bool

	for i := 0; i < len(aln[0].Seq); i++ {
		diff = false
		for j := 1; j < len(aln); j++ {
			if aln[j].Seq[i] != aln[0].Seq[i] {
				diff = true
			}
		}
		if diff {
			S++
		}
	}
	return S
}

//calculateTajimaK is a helper function of Tajima that returns the mean pairwise distance between sequences in an input multiFa alignment block.
func calculateTajimaK(aln []*fasta.Fasta) float64 {
	var distList []int
	for i := 0; i < len(aln); i++ {
		for j := i + 1; j < len(aln); j++ {
			curr := dna.Dist(aln[i].Seq, aln[j].Seq)
			distList = append(distList, curr)
		}
	}

	n := len(distList)
	var sum int
	for i := 0; i < len(distList); i++ {
		sum += distList[i]
	}

	return float64(sum) / float64(n)
}
