package main

import (
	"fmt"
	"github.com/vertgenlab/gonomics/chromInfo"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/numbers"
	"github.com/vertgenlab/gonomics/sam"
	"github.com/vertgenlab/gonomics/vcf"
	"sort"
	"strings"
)

// getVariant finds variants between the experimental and normal samples with p < maxP via fishers exact test.
// returns a vcf record for the input position and a bool that is false if the returned variant should be discarded.
func getVariant(exp, norm []sam.Pile, ref *fasta.Seeker, chroms []chromInfo.ChromInfo, maxP float64, minAf float64, minCoverage int) (vcf.Vcf, bool) {
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

	// pull the base before, as well as the current base in case an anchor is need for an indel
	var seekStart, seekEnd int // Note that the Pos in pile is 1-base
	seekStart = int(bkgd.Pos) - 2
	seekEnd = int(bkgd.Pos)
	if seekStart == -1 {
		seekStart = 0 // TODO get some logic to add the base after a start-of-chromosome deletion for anchor
	}
	refBases, err := fasta.SeekByIndex(ref, bkgd.RefIdx, seekStart, seekEnd) // The Pos in pile is 1-base
	dna.AllToUpper(refBases)                                                 // very important to get the correct index in Pile.Count since lowercase is not allowed in sam/bam
	exception.PanicOnErr(err)
	if len(refBases) == 1 { // only happens if seekStart was negative to begin with //TODO this is lame, fix
		// prepend an N so indexing is as expected
		refBases = append([]dna.Base{dna.N}, refBases...)
	}

	var possibleAlts, passingAlts []string
	var isInsertion []bool                                        // len == len(possibleAlts), true if passingAlts[i] is insertion
	possibleAlts, isInsertion = getPossibleAlts(exp, refBases[1]) // determine alts with min 1 alt read
	altPvalues := make([][]float64, len(possibleAlts))            // first slice is possibleAlts, second slice is experimental samples
	for i := range altPvalues {
		altPvalues[i] = make([]float64, len(exp))
	}
	for i := range possibleAlts {
		for j := range exp {
			if exp[j].RefIdx == -1 {
				altPvalues[i][j] = 1
				continue
			}
			altPvalues[i][j] = fishersExactTest(possibleAlts[i], exp[j], bkgd, hasNorm, minAf, minCoverage, isInsertion[i])
		}
	}

	var passingAltPvalues [][]float64 // subset of altPvalues returned from getPassingAlts
	var passingIsInsertion []bool     // which indexes in passingAlts are insertions
	passingAlts, passingAltPvalues, passingIsInsertion = getPassingAlts(possibleAlts, altPvalues, isInsertion, maxP)
	if len(passingAlts) == 0 {
		return vcf.Vcf{}, false
	}

	// assemble vcf
	var v vcf.Vcf
	v.Samples = make([]vcf.Sample, len(exp)+len(norm))
	v.Chr = chroms[bkgd.RefIdx].Name
	v.Pos = int(bkgd.Pos)
	v.Filter = strings.Join(warnings, ";")
	v.Id = "."
	v.Format = []string{"GT", "DP", "AD", "PV"} // GT = genotype, DP = total reads, AD = reads per allele, PV = p value
	v.Info = "."
	v = addFmtField(v, exp, norm, refBases[1], passingAlts, passingAltPvalues, passingIsInsertion)
	v.Ref = dna.BaseToString(refBases[1]) // an anchor base may be prepended later in the case of deletion
	v.Alt = passingAlts

	var containsDeletion bool
	for i := range v.Alt {
		if v.Alt[i] == "-" {
			containsDeletion = true
		}
		if passingIsInsertion[i] {
			v.Alt[i] = v.Ref + v.Alt[i]
		}
	}

	if containsDeletion {
		// TODO get some logic to add the base after a start-of-chromosome deletion for anchor
		v.Ref = dna.BaseToString(refBases[0]) + v.Ref // prepend anchor base
		for i := range v.Alt {
			if v.Alt[i] == "-" { // replace with previous base
				v.Alt[i] = v.Ref[:1]
			} else { // prepend previous base
				v.Alt[i] = v.Ref[:1] + v.Alt[i]
			}
		}
	}

	return v, true
}

// addFmtField assembles the piles and p values into the format field of the vcf
func addFmtField(v vcf.Vcf, exp, norm []sam.Pile, ref dna.Base, alts []string, passingAltPvalues [][]float64, passingIsInsertion []bool) vcf.Vcf {
	var genotypeAlleles []int16
	var depth int
	var alleleCounts []int
	var pVals []float64

	allSamples := append(exp, norm...)
	for i := range allSamples {
		genotypeAlleles, depth, alleleCounts, pVals = getFormatData(allSamples[i], i, ref, alts, passingAltPvalues, passingIsInsertion)
		if i >= len(exp) { // is normal
			pVals = []float64{-1}
		}
		v.Samples[i].Alleles = genotypeAlleles
		v.Samples[i].Phase = make([]bool, len(genotypeAlleles))
		v.Samples[i].FormatData = []string{"", fmt.Sprintf("%d", depth), sprintAD(alleleCounts), sprintPV(pVals)}
	}
	return v
}

// calcDepth returns the number of reads in the input pile
func calcDepth(s sam.Pile) int {
	var depth int
	for i := range s.Count {
		depth += s.Count[i]
	}
	for _, val := range s.InsCount {
		depth += val
	}
	return depth
}

// getFormatData gathers all information for the FORMAT field of a vcf record.
func getFormatData(s sam.Pile, sIdx int, ref dna.Base, alts []string, passingAltPvalues [][]float64, passingIsInsertion []bool) (genotypeAlleles []int16, depth int, alleleCounts []int, pVals []float64) {
	depth = calcDepth(s)
	pVals = make([]float64, len(alts))

	// add ref to altcounts
	alleleCounts = append(alleleCounts, s.Count[int(ref)])

	for i := range alts {
		if !passingIsInsertion[i] { // not insertion
			alleleCounts = append(alleleCounts, s.Count[int(dna.RuneToBase(rune(alts[i][0])))])
		} else { // is insertion
			alleleCounts = append(alleleCounts, s.InsCount[alts[i]])
		}
		if sIdx < len(passingAltPvalues[i]) { // get p-value only if it is an experimental sample
			pVals[i] = passingAltPvalues[i][sIdx]
		}
	}

	for i := range alleleCounts {
		if alleleCounts[i] > 0 {
			genotypeAlleles = append(genotypeAlleles, int16(i))
		}
	}
	if len(genotypeAlleles) == 1 { // make homozygotes diploid //TODO kind of a hack to assume diploidy, how do others determine ploidy
		genotypeAlleles = append(genotypeAlleles, genotypeAlleles[0])
	}
	return
}

// sprintAD prints the allele depth (AD) field to a string
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

// sprintPV prints the p values (PV) field to a string
func sprintPV(pv []float64) string {
	if len(pv) == 1 && pv[0] == -1 {
		return "."
	}
	var s strings.Builder
	for i := range pv {
		if i > 0 {
			s.WriteByte(',')
		}
		s.WriteString(fmt.Sprintf("%.0g", pv[i]))
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
func getPossibleAlts(exp []sam.Pile, ref dna.Base) ([]string, []bool) {
	if len(exp) == 0 {
		return nil, nil
	}
	var possibleAlts []string
	var isInsertion []bool
	sum := sumPiles(exp)
	for i := range sum.Count { // find possible bases
		if i == int(ref) { // ignore ref base
			continue
		}
		if sum.Count[i] > 0 {
			possibleAlts = append(possibleAlts, string(dna.BaseToRune(dna.Base(i))))
			isInsertion = append(isInsertion, false)
		}
	}

	insertionStartIdx := len(possibleAlts)
	for key, count := range sum.InsCount { // find possible insertions
		if count > 0 {
			possibleAlts = append(possibleAlts, key)
			isInsertion = append(isInsertion, true)
		}
	}

	// sort if multiple insertion sequences are present for stable output
	if len(possibleAlts) > insertionStartIdx + 1 {
		toSort := possibleAlts[insertionStartIdx:]
		sort.Slice(toSort, func(i, j int) bool {
			if len(toSort[i]) == len(toSort[j]) {
				return toSort[i] < toSort[j]
			}
			return len(toSort[i]) < len(toSort[j])
		})
	}

	return possibleAlts, isInsertion
}

// getPassingAlts returns a slice of strings with all alternate alleles with at least 1 samples with a significant p values
// as well as a slice of p values for all samples for each passing alternate allele.
func getPassingAlts(possibleAlts []string, altPvalues [][]float64, isInsertion []bool, maxP float64) (passingAlts []string, passingAltPvalues [][]float64, passingIsInsertion []bool) {
	var pass bool
	passingIsInsertion = make([]bool, 0, len(isInsertion))
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
			passingIsInsertion = append(passingIsInsertion, isInsertion[i])
		}
	}
	return
}

// fishersExactTest returns the p value for the alternate allele present in the exp pile, based on the bkgd pile.
func fishersExactTest(altString string, exp sam.Pile, bkgd sam.Pile, hasNorm bool, minAf float64, minCoverage int, isInsertion bool) float64 {
	// Begin gathering parameters for Fishers Exact Test done in the numbers package
	// test is for the matrix:
	// [a b]
	// [c d]
	// a = Experimental Depth
	// b = Background Depth
	// c = Experimental Alt Allele Count
	// d = Background Alt Allele Count

	var a, b, c, d int

	if !isInsertion { // alt in single base
		alt := int(dna.RuneToBase(rune(altString[0])))
		c = exp.Count[alt]
		d = bkgd.Count[alt]
	} else { // alt in insertion
		c = exp.InsCount[altString]
		d = bkgd.InsCount[altString]
	}

	a = calcDepth(exp) - c
	b = calcDepth(bkgd) - d

	// if no normal is present then bkgd is the sum of all experimental bases,
	// therefore the current experimental base must be subtracted from bkgd totals.
	if !hasNorm {
		b -= a
		d -= c
	}

	// Check exclusion cases to avoid actually doing the test
	var p float64
	switch {
	// If alternate allele is zero then there is no variant and score is 1
	case c == 0:
		p = 1

	// If a = b and c = d then it is testing itself and should return 1
	case a == b && c == d:
		p = 1

	// If the allele frequency of d > c then p is 1
	case float64(c)/float64(c+a) < float64(d)/float64(d+b):
		p = 1

	// If the coverage is less then minCoverage then p is noted as 1 so as to be excluded
	case a+c < minCoverage:
		p = 1

	// If the allele frequency is less than the threshold then p is noted as 1 so as to be excluded
	case float64(c)/float64(c+a) < minAf:
		p = 1

	// If no exclusion conditions are met, then calculate p value
	default:
		if a < c && b < d { // alt is common allele
			p = numbers.FisherExact(a, b, c, d, false)
		} else { // alt is rare uncommon allele
			p = numbers.FisherExact(a, b, c, d, true)
		}
	}
	return p
}
