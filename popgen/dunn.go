package popgen

import (
	//"log".
	"math"

	//"fmt".
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

//below is an implementation of the calculations for Tajima's D. Not currently used in the
//code base, so currently commented out.

/*
//TajimaFromBed calculates Tajima's D from a multiFa alignment at a position set specific by an input bed entry.
//Requires a gapless alignment.
func TajimaFromBedNoGroup(b bed.Bed, aln []fasta.Fasta) float64 {
	bLen := b.ChromEnd - b.ChromStart
	alnPos := fasta.RefPosToAlnPos(aln[0], int(b.ChromStart))
	tmpFa := fasta.CopySubset(aln, alnPos, alnPos+int(bLen))
	tmpFa = fasta.RemoveMissingMult(tmpFa)

	return Tajima(tmpFa)
}

//TajimaFromBed caculates Tajima's D from a multiFa alignment at a position set specified by an input bed entry.
//This version considers only the alignment entries that are represented in an input Group slice. A second return from this function lists members of the input Group slice not found in the multiFa alignment.
//Requires a gapless alignment.
func TajimaFromBed(b bed.Bed, aln []fasta.Fasta, g []*Group) (float64, string) {
	bLen := b.ChromEnd - b.ChromStart
	alnPos := fasta.RefPosToAlnPos(aln[0], int(b.ChromStart))
	tmpFa := fasta.CopySubset(aln, alnPos, alnPos+int(bLen))
	tmp2Fa := fasta.RemoveMissingMult(tmpFa)

	tmp3Fa := FilterMultByGroup(tmp2Fa, g)
	missing := FindMissingGroupMembers(tmp3Fa, g)

	return Tajima(tmpFa), missing
}

//Tajima calculates Tajima's D from an input multiFa alignment block.
func Tajima(aln []fasta.Fasta) float64 {
	k := calculateTajimaK(aln)
	SoverA := calculateTajimaSoverA(aln)
	D := (k - SoverA) / calculateTajimaDenominator(aln)
	return D
}

//calculateTajimaDenominator is a helper function of Tajima that calculates the denominator of the expression.
func calculateTajimaDenominator(aln []fasta.Fasta) float64 {
	n := float64(len(aln))
	b1 := (n + 1) / (3 * (n - 1))
	b2 := 2 * (math.Pow(n, 2) + n + 3) / (9 * n * (n - 1))
	a1 := calculateTajimaA1(aln)
	a2 := calculateTajimaA2(aln)

	c1 := b1 - (1 / a1)
	e1 := c1 / a1
	c2 := b2 - (n+2)/(a1*n) + a2/math.Pow(a1, 2)
	e2 := c2 / (math.Pow(a1, 2) + a2)

	S := float64(calculateS(aln))
	return math.Sqrt(e1*S + e2*S*(S-1))
}

//calculateTajimaA2 is a helper function of calculateDenominator that returns the value of A2, one parameter of Tajima's D.
func calculateTajimaA2(aln []fasta.Fasta) float64 {
	var a2 float64

	for i := 1; i < (len(aln) - 1); i++ {
		a2 += (1.0 / math.Pow(float64(i), 2))
	}

	return a2
}

//calculateTajimaA1 is a helper function of calculateDenominator that returns the value of the A1 parameter.
func calculateTajimaA1(aln []fasta.Fasta) float64 {
	var a1 float64

	for i := 1; i < len(aln)-1; i++ {
		a1 += (1.0 / float64(i))
	}
	return a1
}

//calculateTajimaSoverA is a helper function of Tajima that returns the SoverA approximation of theta.
func calculateTajimaSoverA(aln []fasta.Fasta) float64 {
	S := calculateTajimaS(aln)
	a1 := calculateTajimaA1(aln)
	return (float64(S) / a1)
}

//calculateTajimaS returns the number of segregating sites from a multiFa alignment block.
func calculateTajimaS(aln []fasta.Fasta) int {
	S := 0
	var diff bool

	for i := range aln[0].Seq {
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
func calculateTajimaK(aln []fasta.Fasta) float64 {
	var distList []int
	for i := 0; i < len(aln); i++ {
		for j := i + 1; j < len(aln); j++ {
			curr := dna.Dist(aln[i].Seq, aln[j].Seq)
			distList = append(distList, curr)
		}
	}

	n := len(distList)
	var sum int
	for i := range distList {
		sum += distList[i]
	}

	return float64(sum) / float64(n)
}
*/

//Below are the testing files for the unused Tajima's D implementation.
//I'm not certain of the accuracy of this program, so I would suggest additional tests if anyone wants to pick this back up.
//-Riley
/*
package popgen

import (
	"fmt"
	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/fasta"
	"testing"
)

var seqA []dna.Base = dna.StringToBases("ATAATAAAAAAATAATAAAAAAATAAAAAAAATAAAAAAA")
var seqB []dna.Base = dna.StringToBases("AAAAAAAATAAATAATAAAAAAATAAAAAAAAAAAAAAAA")
var seqC []dna.Base = dna.StringToBases("AAAATAAAAATATAATAAAAAAATATAAAAAAAAAAAAAA")
var seqD []dna.Base = dna.StringToBases("AAAAAAAAAAAATAATAAAAAAATAAATAAATAAAAAAAA")
var seqE []dna.Base = dna.StringToBases("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA")
var seqF []dna.Base = dna.StringToBases("AAAAAAAAAATATAATAAAAAAATAAAAAAAAAAAAAAAA")
var seqG []dna.Base = dna.StringToBases("AAAAAAAAAATATAAAAAAAAAATAAAAAAAAAATTAAAA")
var seqH []dna.Base = dna.StringToBases("AAAAAAAAAATATAATAAAAAAATAAAAAAAAAAATAAAA")
var testFa []fasta.Fasta = []fasta.Fasta{{"eggplant", seqA}, {"raddish", seqB}, {"rhubarb", seqC}, {"asparagus", seqD}, {"broccoli", seqE}, {"tomato", seqF}, {"celery", seqG}, {"carrot", seqH}}
var expected float64 = -1.296575

var b []*bed.Bed = []*bed.Bed{{Chrom: "chr10", ChromStart: 49396820, ChromEnd: 68756350}, {Chrom: "chr10", ChromStart: 75967636, ChromEnd: 76282688}}

func TestTajima(t *testing.T) {
	input := Tajima(testFa)
	//fmt.Printf("%f\n",input)
	if fmt.Sprintf("%f", input) != fmt.Sprintf("%f", expected) {
		t.Errorf("Do not match. Input: %f. Expected: %f.", input, expected)
	}
}

func TestVCFTajima(t *testing.T) {
	input := TajimaGVCFBedSet(b, "testdata/tajimaSet200.vcf")
	expected := -0.998926
	if fmt.Sprintf("%f", input) != fmt.Sprintf("%f", expected) {
		t.Errorf("Do not match. Input: %f. Expected: %f.", input, expected)
	}
}
*/

//Here's the outdata gVCF implementation of Tajima's D as well.
/*
package popgen

import (
	"fmt"
	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/interval"
	"github.com/vertgenlab/gonomics/vcf"
	"log"
	//DEBUG: "fmt"
	"math"
	"strings"
)

This file contains functions for calculating Tajima's D from a gVCF. More information on Tajima's D can be found at https://en.wikipedia.org/wiki/Tajima%27s_D

//IndividualAllele is an internal struct for the Tajima's D calculation from gVCFs. It contains the "column" information of a gVCF block
//In other words, an IndividualAllele represents the state of each segregating site for one allele in a gVCF block.
type IndividualAllele struct {
	sites [][]dna.Base
}

//IndividualDist is a helper function for Tajima's D calculations that calculates the distance (the number of segregating sites that are different by state) between two IndividualAllele structs.
func IndividualDist(a IndividualAllele, b IndividualAllele) int {
	var answer int = 0
	if len(a.sites) != len(b.sites) {
		log.Fatalf("IndividualAllele sites lists must be the same length to calculate distance.")
	}
	for i := 0; i < len(a.sites); i++ {
		if len(a.sites[i]) != len(b.sites[i]) {
			answer++
		} else if dna.Dist(a.sites[i], b.sites[i]) > 0 { //dist requires sequences of equal length, so the first if catches this case.
			answer++
		}
	}
	return answer
}

//GVCFCalculateA2 returns the A2 sum used in the calculation of Tajima's D. More information for this constant can be found at the wikipedia page.
func GVCFCalculateA2(n int) float64 {
	var answer float64

	for i := 1; i < n; i++ {
		answer += 1.0 / math.Pow(float64(i), 2)
	}
	return answer
}

//GVCFCalculateA1 returns the A1 sum used in the calculation of Tajima's D. (see wiki)
func GVCFCalculateA1(n int) float64 {
	var answer float64

	for i := 1; i < n; i++ {
		answer += 1.0 / float64(i)
	}
	return answer
}

//CalculatePiTajima calculates Pi, the mean pairwise distance between the members of an individual allele matrix.
func CalculatePiTajima(all []*IndividualAllele) float64 {
	var sumDiffs, numComparisons int
	for i := 0; i < len(all); i++ {
		for j := i + 1; j < len(all); j++ { //for all pairs of alleles
			sumDiffs += IndividualDist(*all[i], *all[j])
			numComparisons++
		}
	}
	return float64(sumDiffs) / float64(numComparisons)
}

//CountSegSites returns the number of segregating sites in a slice of IndividualAllele structs.
func CountSegSites(all []*IndividualAllele) int {
	var segSites int = 0
	var found bool = false

	if len(all) < 2 {
		log.Fatalf("Need more than two alleles to count segregating sites\n")
	}

	for s := range all[0].sites {
		found = false
		for i := 0; i < len(all) && found == false; i++ {
			for j := i + 1; j < len(all) && found == false; j++ {
				if len(all[i].sites[s]) != len(all[j].sites[s]) || dna.Dist(all[i].sites[s], all[j].sites[s]) > 0 {
					segSites++
					found = true
				}
			}
		}
	}
	return segSites
}

//TajimaGVCF is the main function that takes an individual allele matrix and calculates Tajima's D. Limited to SNPs.
func TajimaGVCF(all []*IndividualAllele) float64 {
	pi := CalculatePiTajima(all)
	n := len(all)
	S := CountSegSites(all)
	a1 := GVCFCalculateA1(n)
	a2 := GVCFCalculateA2(n)
	b1 := (float64(n + 1)) / (3.0 * float64(n-1))
	b2 := 2.0 * (math.Pow(float64(n), 2) + float64(n) + 3.0) / (9.0 * float64(n) * float64(n-1))
	c1 := b1 - (1.0 / a1)
	e1 := c1 / a1
	c2 := b2 - float64(n+2)/(a1*float64(n)) + a2/math.Pow(a1, 2)
	e2 := c2 / (math.Pow(a1, 2) + a2)
	return (pi - (float64(S) / a1)) / math.Sqrt(e1*float64(S)+e2*float64(S)*float64(S-1))
}

//TajimaGVCFBedSet is designed to calculate Tajima's D from a
//gVCF file for variants falling in a particular set of bed regions. Returns one value for the whole set. Limited to SNPs.
func TajimaGVCFBedSet(b []*bed.Bed, VcfFile string) float64 {
	alpha, _ := vcf.GoReadToChan(VcfFile)
	var all []*IndividualAllele
	var j, k int
	var firstTime bool = true
	var intervals []interval.Interval
	intervals = make([]interval.Interval, len(b))
	for i := range b {
		intervals[i] = b[i]
	}
	tree := interval.BuildTree(intervals)
	for i := range alpha {
		if len(interval.Query(tree, &i, "any")) > 0 { //in other words, if the variant overlaps any of the beds
			if !strings.ContainsAny(i.Alt[0], "<>") { //We convert the alt and ref to []DNA.base, so structural variants with <CN0> notation will fail to convert. This check allows us to ignore these cases.
				//DEBUG: fmt.Printf("Len: %v.\n", len(g.Genotypes))
				if firstTime {
					firstTime = false
					all = make([]*IndividualAllele, len(i.Samples)*2) //makes the individual list, with one entry for each allele

					for k = range all { //in this loop we initialize all the sites lists
						currSites := make([][]dna.Base, 0)
						all[k] = &IndividualAllele{sites: currSites}
					}
				}
				for j = range i.Samples {
					if len(i.Samples[j].Alleles) != 2 || i.Samples[j].Alleles[0] == -1 || i.Samples[j].Alleles[1] == -1 { //check that data exists for both alleles
						fmt.Println(i)
						log.Fatalf("Tajima's D on VCFs requires complete alignment blocks.")
					} else {
						if i.Samples[j].Alleles[0] == 0 {
							all[2*j].sites = append(all[2*j].sites, dna.StringToBases(i.Ref))
						} else {
							all[2*j].sites = append(all[2*j].sites, dna.StringToBases(i.Alt[0])) //TODO: (riley) currently only handles biallelic positions.
						}
						if i.Samples[j].Alleles[1] == 0 {
							all[2*j+1].sites = append(all[2*j+1].sites, dna.StringToBases(i.Ref))
						} else {
							all[2*j+1].sites = append(all[2*j+1].sites, dna.StringToBases(i.Alt[0]))
						}
					}
				}
			}
		}
	}
	return TajimaGVCF(all)
}

*/
