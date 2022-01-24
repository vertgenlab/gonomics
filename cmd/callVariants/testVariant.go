package main

import (
	"fmt"
	"github.com/vertgenlab/gonomics/chromInfo"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/sam"
	"github.com/vertgenlab/gonomics/vcf"
	"strings"
)

func getVariant(exp, norm []sam.Pile, ref *fasta.Seeker, chroms []chromInfo.ChromInfo, maxP float64) (vcf.Vcf, bool) {
	var hasExp, hasNorm bool
	var warnings []string
	for i := range exp {
		if exp[i].RefIdx != -1 { // RefIdx == -1 is mark for no data
			hasExp = true
		}
	}
	if !hasExp { // return empty struct if there are no experimental samples for the input site
		return vcf.Vcf{}, false
	}
	for i := range norm {
		if norm[i].RefIdx != -1 {
			hasNorm = true
		}
	}
	if !hasNorm {
		warnings = append(warnings, "NO_NORMAL")
	}

	// get a background allele count for fishers exact testing
	var bkgd sam.Pile
	if hasNorm { // if normals are present, then use the summed normals for comparison
		bkgd = sumPiles(norm)
	} else { // if no normals are present then use experimental samples for comparison
		bkgd = sumPiles(exp)
	}

	refBases, err := fasta.SeekByIndex(ref, bkgd.RefIdx, int(bkgd.Pos), int(bkgd.Pos+1)) // pull the base before and after //TODO double check this
	exception.PanicOnErr(err)

	var possibleAlts, passingAlts []string
	possibleAlts = getPossibleAlts(exp, refBases[0]) // determine alts with min 1 alt read
	altPvalues := make([][]float64, len(possibleAlts)) // first slice is possibleAlts, second slice is experimental samples
	for i := range altPvalues {
		altPvalues[i] = make([]float64, len(exp))
	}
	for i := range possibleAlts {
		for j := range exp {
			if exp[j].RefIdx == -1 {
				altPvalues[i][j] = 1
				continue
			}
			altPvalues[i][j] = fishersExactTest(possibleAlts[i], exp[j], bkgd, hasNorm)
		}
	}

	var passingAltPvalues [][]float64 // subset of altPvalues returned from getPassingAlts
	passingAlts, passingAltPvalues = getPassingAlts(possibleAlts, altPvalues, maxP)
	if len(passingAlts) == 0 {
		return vcf.Vcf{}, false
	}

	var v vcf.Vcf
	v.Samples = make([]vcf.Sample, len(exp) + len(norm))
	v.Chr = chroms[bkgd.RefIdx].Name
	v.Pos = int(bkgd.Pos) + 1
	v.Filter = strings.Join(warnings, ";")
	v.Id = "."
	v.Format = []string{"GT", "DP", "AD", "PV"} // GT = genotype, DP = total reads, AD = reads per allele, PV = p value
	v.Info = "."
	v.Alt = passingAlts
	v.Ref = dna.BaseToString(dna.ToUpper(refBases[0])) //TODO check for indels

	return makeVcf(v, exp, norm, bkgd, hasNorm, refBases[0], passingAlts, passingAltPvalues), true
}

func makeVcf(v vcf.Vcf, exp, norm []sam.Pile, bkgd sam.Pile, hasNorm bool, ref dna.Base, alts []string, altPvalues [][]float64) vcf.Vcf {
	var genotypeAlleles []int16
	var depth int
	var alleleCounts []int
	var pVal float64

	allSamples := append(exp, norm...)
	for i := range allSamples {
		genotypeAlleles, depth, alleleCounts, pVal = getFormatData(allSamples[i], bkgd, i, hasNorm, ref, alts, altPvalues)
		if i >= len(exp) { // is normal
			pVal = 1
		}
		v.Samples[i].Alleles = genotypeAlleles
		v.Samples[i].Phase = make([]bool, len(genotypeAlleles))
		v.Samples[i].FormatData = []string{"", fmt.Sprintf("%d", depth), sprintAD(alleleCounts), fmt.Sprintf("%.0E", pVal)}
	}
	return v
}

//TODO provides dummy values atm
func getFormatData(s, bkgd sam.Pile, sIdx int, hasNorm bool, ref dna.Base, alts []string, altPvalues [][]float64) (genotypeAlleles []int16, depth int, alleleCounts []int, pVal float64) {
	for i := range s.Count {
		depth += s.Count[i]
	}
	for _, val := range s.InsCount {
		depth += val
	}

	// add ref to altcounts
	alleleCounts = append(alleleCounts, s.Count[int(ref)])

	for i := range alts {
		if len(alts[i]) == 1 { // not insertion
			alleleCounts = append(alleleCounts, s.Count[int(dna.RuneToBase(rune(alts[i][0])))])
		} else { // is insertion
			alleleCounts = append(alleleCounts, s.InsCount[alts[i]])
		}
	}

	for i := range alleleCounts {
		if alleleCounts[i] > 0 {
			genotypeAlleles = append(genotypeAlleles, int16(i))
		}
	}
	if len(genotypeAlleles) == 1 { // make homozygotes  diploid //TODO kind of a hack to assume diploidy, how do others determine ploidy
		genotypeAlleles = append(genotypeAlleles, genotypeAlleles[0])
	}

	pVal = 0.01
	return
}

func sprintAD(ad []int) string {
	var s strings.Builder
	for i := range ad {
		if i > 0 {
			s.WriteByte(',')
		}
		s.WriteString(fmt.Sprintf("%d", ad[i]))
	}
	return s.String()
}

// sumPiles outputs a Pile with counts equal to the sum of all per-base counts of Piles in p.
func sumPiles(p []sam.Pile) sam.Pile {
	var answer sam.Pile
	answer.InsCount = make(map[string]int)
	answer.RefIdx = -1

	for i := range p {
		if answer.RefIdx == -1 && p[i].RefIdx != -1 {
			answer.RefIdx = p[i].RefIdx
			answer.Pos = p[i].Pos
		}
		for j := range answer.Count {
			answer.Count[j] += p[i].Count[j]
		}
		for key, val := range p[i].InsCount {
			answer.InsCount[key] += val
		}
	}
	return answer
}

// getPossibleAlts returns a slice of strings containing every non-reference allele with at least 1 supporting read.
func getPossibleAlts(exp []sam.Pile, ref dna.Base) []string {
	if len(exp) == 0 {
		return nil
	}
	var possibleAlts []string
	sum := sumPiles(exp)
	for i := range sum.Count { // find possible bases
		if i == int(ref) { // ignore ref base
			continue
		}
		if sum.Count[i] > 0 {
			possibleAlts = append(possibleAlts, string(dna.BaseToRune(dna.Base(i))))
		}
	}

	for key, count := range sum.InsCount { // find possible insertions
		if count > 0 {
			possibleAlts = append(possibleAlts, key)
		}
	}
	return possibleAlts
}

// getPassingAlts returns a slice of strings with all alternate alleles with at least 1 samples with a significant p values
// as well as a slice of p values for all samples for each passing alternate allele.
func getPassingAlts(possibleAlts []string, altPvalues [][]float64, maxP float64) (passingAlts []string, passingAltPvalues [][]float64) {
	var pass bool
	for i := range altPvalues { // first slice is possibleAlts, second slice is experimental samples
		pass = false
		for j := range altPvalues[i] {
			if altPvalues[i][j] < maxP {
				pass = true
				break
			}
		}
		if pass {
			passingAlts = append(passingAlts, possibleAlts[i])
			passingAltPvalues = append(passingAltPvalues, altPvalues[i])
		}
	}
	return
}

func fishersExactTest(alt string, exp sam.Pile, bkgd sam.Pile, hasNorm bool) float64 {
	return 0.001 //TODO
}