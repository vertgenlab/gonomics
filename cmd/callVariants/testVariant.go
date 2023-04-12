package main

import (
	"fmt"
	"sort"
	"strconv"
	"strings"

	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/numbers"
	"github.com/vertgenlab/gonomics/sam"
	"github.com/vertgenlab/gonomics/vcf"
)

type variantType int

const (
	singleNucleotide = iota
	insertion
	deletion
)

// getVariant finds variants between the experimental and normal samples with p < maxP via fishers exact test.
// returns a vcf record for the input position and a bool that is false if the returned variant should be discarded.
func getVariant(exp, norm []sam.Pile, samHeader sam.Header, ref *fasta.Seeker, maxP, minAf, maxAf, maxStrandBias float64, minCoverage, minAltReads int) (vcf.Vcf, bool) {
	var warnings []string
	var bkgd sam.Pile // background allele count for fishers exact testing
	var hasNorm bool

	if !dataPresent(exp) { // return empty struct if there are no experimental samples for the input site
		return vcf.Vcf{}, false
	}
	if !dataPresent(norm) {
		bkgd = sumPiles(exp) // if no normals are present then use experimental samples for comparison
		warnings = append(warnings, "NO_NORMAL")
	} else {
		bkgd = sumPiles(norm) // if normals are present, then use the summed normals for comparison
		hasNorm = true
	}

	// pull the base before, as well as the current base in case an anchor is need for an indel
	// len(refBases) always == 2
	chrName := samHeader.Chroms[bkgd.RefIdx].Name
	refBases := getRef(int(bkgd.Pos)-1, int(bkgd.Pos), chrName, ref)

	var alts []string
	var altPvalues [][]float64
	var altVarTypes []variantType
	alts, altPvalues, altVarTypes = getAlts(exp, bkgd, refBases, hasNorm, minAf, maxAf, maxStrandBias, minCoverage, minAltReads, maxP)

	if len(alts) == 0 {
		return vcf.Vcf{}, false
	}

	return makeVcf(exp, norm, bkgd, chrName, warnings, refBases, alts, altPvalues, altVarTypes, ref), true
}

// addFmtField assembles the piles and p values into the format field of the vcf
func addFmtField(v vcf.Vcf, exp, norm []sam.Pile, ref dna.Base, alts []string, passingAltPvalues [][]float64, passingVarType []variantType) vcf.Vcf {
	var genotypeAlleles []int16
	var depth int
	var alleleCounts []int
	var pVals []float64

	allSamples := append(exp, norm...)
	for i := range allSamples {
		genotypeAlleles, depth, alleleCounts, pVals = getFormatData(allSamples[i], i, ref, alts, passingAltPvalues, passingVarType)
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
	for i := range s.CountF {
		depth += s.CountF[i] + s.CountR[i]
	}
	for _, val := range s.InsCountF {
		depth += val
	}
	for _, val := range s.InsCountR {
		depth += val
	}
	// Note that DelCount does NOT count towards depth
	return depth
}

// getFormatData gathers all information for the FORMAT field of a vcf record.
func getFormatData(s sam.Pile, sIdx int, ref dna.Base, alts []string, passingAltPvalues [][]float64, passingVarType []variantType) (genotypeAlleles []int16, depth int, alleleCounts []int, pVals []float64) {
	depth = calcDepth(s)
	pVals = make([]float64, len(alts))
	var err error
	var delInt int

	// add ref to altcounts
	alleleCounts = append(alleleCounts, s.CountF[int(ref)]+s.CountR[int(ref)])
	var b dna.Base
	for i := range alts {
		switch passingVarType[i] {
		case singleNucleotide:
			b, err = dna.RuneToBase(rune(alts[i][0]))
			exception.PanicOnErr(err)
			alleleCounts = append(alleleCounts, s.CountF[int(b)]+s.CountR[int(b)])

		case insertion:
			alleleCounts = append(alleleCounts, s.InsCountF[alts[i]]+s.InsCountR[alts[i]])

		case deletion:
			delInt, err = strconv.Atoi(alts[i]) // this is wicked lame
			exception.PanicOnErr(err)
			alleleCounts = append(alleleCounts, s.DelCountF[delInt]+s.DelCountR[delInt])
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
	if len(genotypeAlleles) == 1 { // make homozygotes diploid //TODO kind of a hack to assume diploidy, how do others determine ploidy?
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
	answer.InsCountF = make(map[string]int)
	answer.InsCountR = make(map[string]int)
	answer.DelCountF = make(map[int]int)
	answer.DelCountR = make(map[int]int)
	answer.RefIdx = -1

	for i := range p {
		if answer.RefIdx == -1 && p[i].RefIdx != -1 {
			answer.RefIdx = p[i].RefIdx
			answer.Pos = p[i].Pos
		}
		for j := range answer.CountF {
			answer.CountF[j] += p[i].CountF[j]
			answer.CountR[j] += p[i].CountR[j]
		}
		for key, val := range p[i].InsCountF {
			answer.InsCountF[key] += val
		}
		for key, val := range p[i].InsCountR {
			answer.InsCountR[key] += val
		}
		for key, val := range p[i].DelCountF {
			answer.DelCountF[key] += val
		}
		for key, val := range p[i].DelCountR {
			answer.DelCountR[key] += val
		}
	}
	return answer
}

// getAlts returns a slice of all qualifying alternate alleles, their p values, and whether they are an insertion.
func getAlts(exp []sam.Pile, bkgd sam.Pile, refBases []dna.Base, hasNorm bool, minAf, maxAf, maxStrandBias float64, minCoverage, minAltReads int, maxP float64) (passingAlts []string, passingAltPvalues [][]float64, variantTypes []variantType) {
	var possibleAlts []string
	possibleAlts, variantTypes = getPossibleAlts(exp, refBases[1]) // determine alts with min 1 alt read
	altPvalues := make([][]float64, len(possibleAlts))             // first slice is possibleAlts, second slice is experimental samples
	for i := range altPvalues {
		altPvalues[i] = make([]float64, len(exp))
	}
	for i := range possibleAlts {
		for j := range exp {
			if exp[j].RefIdx == -1 {
				altPvalues[i][j] = 1
				continue
			}
			altPvalues[i][j] = fishersExactTest(possibleAlts[i], exp[j], bkgd, hasNorm, minAf, maxAf, maxStrandBias, minCoverage, minAltReads, variantTypes[i])
		}
	}

	return getPassingAlts(possibleAlts, altPvalues, variantTypes, maxP)
}

// getPossibleAlts returns a slice of strings containing every non-reference allele with at least 1 supporting read.
func getPossibleAlts(exp []sam.Pile, ref dna.Base) ([]string, []variantType) {
	if len(exp) == 0 {
		return nil, nil
	}
	var possibleAlts []string
	var variantTypes []variantType
	sum := sumPiles(exp)
	for i := range sum.CountF { // find possible bases
		if i == int(ref) || i == int(dna.Gap) { // ignore ref base and dna.Gap
			continue
		}
		if sum.CountF[i] > 0 || sum.CountR[i] > 0 {
			possibleAlts = append(possibleAlts, string(dna.BaseToRune(dna.Base(i))))
			variantTypes = append(variantTypes, singleNucleotide)
		}
	}

	for key, count := range sum.DelCountF { // find possible deletions
		if count > 0 {
			possibleAlts = append(possibleAlts, fmt.Sprintf("%d", key)) // lame
			variantTypes = append(variantTypes, deletion)
		}
	}

	var presentInForward bool
	for key, count := range sum.DelCountR {
		_, presentInForward = sum.DelCountF[key]
		if !presentInForward && count > 0 {
			possibleAlts = append(possibleAlts, fmt.Sprintf("%d", key)) // lame
			variantTypes = append(variantTypes, deletion)
		}
	}

	insertionStartIdx := len(possibleAlts)
	for key, count := range sum.InsCountF { // find possible insertions
		if count > 0 {
			possibleAlts = append(possibleAlts, key)
			variantTypes = append(variantTypes, insertion)
		}
	}

	for key, count := range sum.InsCountR {
		_, presentInForward = sum.InsCountF[key]
		if !presentInForward && count > 0 {
			possibleAlts = append(possibleAlts, key)
			variantTypes = append(variantTypes, insertion)
		}
	}

	// sort if multiple insertion sequences are present for stable output
	if len(possibleAlts) > insertionStartIdx+1 {
		toSort := possibleAlts[insertionStartIdx:]
		sort.Slice(toSort, func(i, j int) bool {
			if len(toSort[i]) == len(toSort[j]) {
				return toSort[i] < toSort[j]
			}
			return len(toSort[i]) < len(toSort[j])
		})
	}

	return possibleAlts, variantTypes
}

// getPassingAlts returns a slice of strings with all alternate alleles with at least 1 samples with a significant p values
// as well as a slice of p values for all samples for each passing alternate allele.
func getPassingAlts(possibleAlts []string, altPvalues [][]float64, variantTypes []variantType, maxP float64) (passingAlts []string, passingAltPvalues [][]float64, passingVariantTypes []variantType) {
	var pass bool
	passingVariantTypes = make([]variantType, 0, len(variantTypes))
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
			passingVariantTypes = append(passingVariantTypes, variantTypes[i])
		}
	}
	return
}

// fishersExactTest returns the p value for the alternate allele present in the exp pile, based on the bkgd pile.
func fishersExactTest(altString string, exp sam.Pile, bkgd sam.Pile, hasNorm bool, minAf, maxAf, maxStrandBias float64, minCoverage, minAltReads int, varType variantType) float64 {
	// Begin gathering parameters for Fishers Exact Test done in the numbers package
	// test is for the matrix:
	// [a b]
	// [c d]
	// a = Experimental Depth
	// b = Background Depth
	// c = Experimental Alt Allele Count
	// d = Background Alt Allele Count

	var a, b, c, d int
	var base dna.Base
	var err error
	var fwdStrandBias float64

	switch varType {
	case singleNucleotide: // alt is single base
		base, err = dna.RuneToBase(rune(altString[0]))
		exception.PanicOnErr(err)
		alt := int(base)
		c = exp.CountF[alt] + exp.CountR[alt]
		d = bkgd.CountF[alt] + bkgd.CountR[alt]
		fwdStrandBias = float64(exp.CountF[alt]) / float64(c)

	case insertion:
		c = exp.InsCountF[altString] + exp.InsCountR[altString]
		d = bkgd.InsCountF[altString] + bkgd.InsCountR[altString]
		fwdStrandBias = float64(exp.InsCountF[altString]) / float64(c)

	case deletion:
		delInt, err := strconv.Atoi(altString) // this is wicked lame
		exception.PanicOnErr(err)
		c = exp.DelCountF[delInt] + exp.DelCountR[delInt]
		d = bkgd.DelCountF[delInt] + bkgd.DelCountR[delInt]
		fwdStrandBias = float64(exp.DelCountF[delInt]) / float64(c)
	}

	if fwdStrandBias > maxStrandBias || fwdStrandBias < 1-maxStrandBias {
		return 1
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

	// If alternate allele is less than minimum reads then p is 1
	case c < minAltReads:
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

	// If the allele frequency is outside the threshold then p is noted as 1 so as to be excluded
	case float64(c)/float64(c+a) < minAf || float64(c)/float64(c+a) > maxAf:
		p = 1

	// If no exclusion conditions are met, then calculate p value
	default:
		p = numbers.FisherExact(a, b, c, d, true)
	}
	return p
}

// dataPresent returns true if any of the input []Pile contain a non-nil value
func dataPresent(p []sam.Pile) bool {
	for i := range p {
		if p[i].RefIdx != -1 { // RefIdx == -1 is mark for no data
			return true
		}
	}
	return false
}

// getRef pulls the requested base and the base before the requested base to use as an anchor for deletions.
// Note that the input pos is zero-based, left-closed, right-open.
func getRef(start, end int, chrName string, ref *fasta.Seeker) []dna.Base {
	var seekStart, seekEnd int
	seekStart = start - 1
	seekEnd = end
	if seekStart == -1 {
		seekStart = 0 // TODO get some logic to add the base after a start-of-chromosome deletion for anchor
	}

	refBases, err := fasta.SeekByName(ref, chrName, seekStart, seekEnd) // The Pos in pile is 1-base
	dna.AllToUpper(refBases)                                            // very important to get the correct index in Pile.Count since lowercase is not allowed in sam/bam
	exception.PanicOnErr(err)

	if len(refBases) == 1 { // only happens if seekStart was negative to begin with
		// prepend an N so indexing is as expected
		refBases = append([]dna.Base{dna.N}, refBases...)
	}
	return refBases
}

// makeVcf asembles a vcf from the input values.
func makeVcf(exp, norm []sam.Pile, bkgd sam.Pile, chrName string, warnings []string, refBases []dna.Base, passingAlts []string, passingAltPvalues [][]float64, passingVarTypes []variantType, ref *fasta.Seeker) vcf.Vcf {
	// assemble vcf
	var v vcf.Vcf
	v.Samples = make([]vcf.Sample, len(exp)+len(norm))
	v.Chr = chrName
	v.Pos = int(bkgd.Pos)
	v.Filter = strings.Join(warnings, ";")
	v.Id = "."
	v.Format = []string{"GT", "DP", "AD", "PV"} // GT = genotype, DP = total reads, AD = reads per allele, PV = p value
	v.Info = "."
	v = addFmtField(v, exp, norm, refBases[1], passingAlts, passingAltPvalues, passingVarTypes)
	v.Ref = dna.BaseToString(refBases[1]) // an anchor base may be prepended later in the case of deletion
	v.Alt = passingAlts                   // will transform deletions later

	deletionIndexes := make([]int, 0, len(passingAlts))
	for i := range passingVarTypes {
		if passingVarTypes[i] == deletion {
			deletionIndexes = append(deletionIndexes, i)
		}
	}

	v = adjustAlts(v, deletionIndexes, passingVarTypes, ref)
	return v
}

// adjustAlts alters the ref and alt fields depending on the variant type
func adjustAlts(v vcf.Vcf, deletionIndexes []int, varTypes []variantType, ref *fasta.Seeker) vcf.Vcf {
	var longestDeletion int
	delLens := make([]int, len(deletionIndexes))
	var err error
	for i := range deletionIndexes {
		delLens[i], err = strconv.Atoi(v.Alt[deletionIndexes[i]])
		exception.PanicOnErr(err)
		if delLens[i] > longestDeletion {
			longestDeletion = delLens[i]
		}
	}

	var hasAnchor bool
	if len(deletionIndexes) > 0 {
		v.Pos-- // since anchor base must be added
		refBases := getRef(v.Pos, v.Pos+longestDeletion, v.Chr, ref)
		v.Ref = dna.BasesToString(refBases)
		hasAnchor = true
	}

	var delLenIdx int
	s := new(strings.Builder)

	for i := range v.Alt {
		switch varTypes[i] {
		case singleNucleotide:
			v.Alt[i] = getSnvAltString(s, v.Ref, v.Alt[i], hasAnchor)

		case insertion:
			v.Alt[i] = getInsAltString(s, v.Ref, v.Alt[i], hasAnchor)

		case deletion:
			v.Alt[i] = getDelAltString(s, v.Ref, delLens[delLenIdx])
			delLenIdx++
		}
	}

	return v
}

func getSnvAltString(s *strings.Builder, ref, alt string, hasAnchor bool) string {
	s.Reset()
	if hasAnchor {
		s.WriteByte(ref[0]) // anchor if present
	}
	s.WriteByte(alt[0]) // alt base
	if len(ref) > 2 {
		s.WriteString(ref[2:]) // any remaining bases
	}
	return s.String()
}

func getInsAltString(s *strings.Builder, ref, insSeq string, hasAnchor bool) string {
	s.Reset()
	if !hasAnchor { // if no anchor, then no deletion is present and simple return
		return ref + insSeq
	}

	// deletion is present, must work around
	s.WriteString(ref[:2] + insSeq) // anchor + refbase + inserted sequence
	if len(ref) > 2 {
		s.WriteString(ref[2:]) // remaining bases
	}
	return s.String()
}

func getDelAltString(s *strings.Builder, ref string, delLen int) string {
	s.Reset()
	// deletions will always have an anchor
	s.WriteByte(ref[0])      // anchor base
	if len(ref) > delLen+1 { // this happens if there is a longer deletion present at the same position
		s.WriteString(ref[delLen+1:]) // the +1 is for the base after the delLen and after the anchor base
	}
	return s.String()
}
