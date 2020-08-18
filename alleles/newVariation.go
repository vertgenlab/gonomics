package alleles

import (
	"fmt"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/numbers"
	"github.com/vertgenlab/gonomics/vcf"
	"path"
	"path/filepath"
)

// FindNewVariation calls variants from a set of samples and normals that DO NOT already exist in the graph structure.
func FindNewVariation(alleleStream <-chan []*Allele, normalIDs map[string]bool, afThreshold float64, sigThreshold float64) <-chan *vcf.Vcf {
	answer := make(chan *vcf.Vcf)

	// To make the normal IDs option if the value is nil, just initialize empty map so all samples
	// are treated as experimental. The call would look like FindNewVariation(alleleStream, nil, afThreshold ...)
	if normalIDs == nil {
		normalIDs = make(map[string]bool)
	}

	go scoreAlleles(answer, alleleStream, normalIDs, afThreshold, sigThreshold)
	return answer
}

// scoreAlleles is Designed to be run as a goroutine that accepts alleles from the alleleStream channel,
// computes the p value and makes a VCF, then sends the vcf record on the answer channel
func scoreAlleles(answer chan<- *vcf.Vcf, alleleStream <-chan []*Allele, normalIDs map[string]bool, afThreshold float64, sigThreshold float64) {

	for alleles := range alleleStream {
		if len(alleles) == 1 {
			continue
		} // Can only evaluate when at least two samples are present
		bkgd, normalPresent := calcBackground(alleles, normalIDs)

		var a, b, c, d int
		var p float64
		var warnings []string

		for i := 0; i < len(alleles); i++ {
			if !normalIDs[alleles[i].Sample] { // if sample is listed as experimental
				//TODO: There is room for improvement here. the b and d values values don't need to be determined
				// for each individual allele and each individual base. We just need to make adjustments on the fly
				// depending on whether or not normals are present and therefore we need to do b -= a and d -= c

				warnings = nil
				if !normalPresent {
					warnings = append(warnings, "NO_NORMAL")
				}

				// Test for A
				a, b, c, d = getFishersInput(alleles[i], bkgd, normalPresent, dna.A, 0)
				p = alleleFishersExact(a, b, c, d, afThreshold)
				if p <= sigThreshold {
					if passesStrandBias(alleles[i].Count.BaseAF, alleles[i].Count.BaseAR) {
						answer <- alleleToVcf(alleles[i], p, dna.A, 0, warnings)
					} else {
						answer <- alleleToVcf(alleles[i], p, dna.A, 0, append(warnings, "STRAND_BIAS"))
					}
				}

				// Test for C
				a, b, c, d = getFishersInput(alleles[i], bkgd, normalPresent, dna.C, 0)
				p = alleleFishersExact(a, b, c, d, afThreshold)
				if p <= sigThreshold {
					if passesStrandBias(alleles[i].Count.BaseCF, alleles[i].Count.BaseCR) {
						answer <- alleleToVcf(alleles[i], p, dna.C, 0, warnings)
					} else {
						answer <- alleleToVcf(alleles[i], p, dna.C, 0, append(warnings, "STRAND_BIAS"))
					}
				}

				// Test for G
				a, b, c, d = getFishersInput(alleles[i], bkgd, normalPresent, dna.G, 0)
				p = alleleFishersExact(a, b, c, d, afThreshold)
				if p <= sigThreshold {
					if passesStrandBias(alleles[i].Count.BaseGF, alleles[i].Count.BaseGR) {
						answer <- alleleToVcf(alleles[i], p, dna.G, 0, warnings)
					} else {
						answer <- alleleToVcf(alleles[i], p, dna.G, 0, append(warnings, "STRAND_BIAS"))
					}
				}

				// Test for T
				a, b, c, d = getFishersInput(alleles[i], bkgd, normalPresent, dna.T, 0)
				p = alleleFishersExact(a, b, c, d, afThreshold)
				if p <= sigThreshold {
					if passesStrandBias(alleles[i].Count.BaseTF, alleles[i].Count.BaseTR) {
						answer <- alleleToVcf(alleles[i], p, dna.T, 0, warnings)
					} else {
						answer <- alleleToVcf(alleles[i], p, dna.T, 0, append(warnings, "STRAND_BIAS"))
					}
				}

				// Test for Indels
				for k := 0; k < len(alleles[i].Count.Indel); k++ {
					a, b, c, d = getFishersInput(alleles[i], bkgd, normalPresent, dna.Gap, k)
					p = alleleFishersExact(a, b, c, d, afThreshold)
					if p <= sigThreshold {
						if passesStrandBias(alleles[i].Count.Indel[k].CountF, alleles[i].Count.Indel[k].CountR) {
							answer <- alleleToVcf(alleles[i], p, dna.Gap, k, warnings)
						} else {
							answer <- alleleToVcf(alleles[i], p, dna.Gap, k, append(warnings, "STRAND_BIAS"))
						}
					}
				}
			}
		}
	}
	close(answer)
}

// calcBackground adds the count values for all allele structs designated as normal and returns a single background struct
func calcBackground(alleles []*Allele, normalIDs map[string]bool) (*Allele, bool) {
	var answer *Allele = &Allele{Count: &AlleleCount{alleles[0].Count.Ref, 0, 0, 0, 0, 0, 0, 0, 0, 0, make([]Indel, 0)}}
	var normalPresent bool

	for i := 0; i < len(alleles); i++ {
		if normalIDs[alleles[i].Sample] {
			addAlleles(answer, alleles[i])
			normalPresent = true
		}
	}

	if normalPresent {
		return answer, normalPresent
	} else {
		for k := 0; k < len(alleles); k++ {
			addAlleles(answer, alleles[k])
		}
	}
	return answer, normalPresent
}

// addAlleles adds the count values in two allele structs
func addAlleles(sum *Allele, toBeAdded *Allele) {
	sum.Count.BaseAF += toBeAdded.Count.BaseAF
	sum.Count.BaseAR += toBeAdded.Count.BaseAR
	sum.Count.BaseCF += toBeAdded.Count.BaseCF
	sum.Count.BaseCR += toBeAdded.Count.BaseCR
	sum.Count.BaseGF += toBeAdded.Count.BaseGF
	sum.Count.BaseGR += toBeAdded.Count.BaseGR
	sum.Count.BaseTF += toBeAdded.Count.BaseTF
	sum.Count.BaseTR += toBeAdded.Count.BaseTR

	for j := 0; j < len(toBeAdded.Count.Indel); j++ {
		indel := findMatchingIndel(&toBeAdded.Count.Indel[j], sum.Count.Indel)
		if indel != nil {
			indel.CountF += toBeAdded.Count.Indel[j].CountF
			indel.CountR += toBeAdded.Count.Indel[j].CountR
		} else {
			sum.Count.Indel = append(sum.Count.Indel, toBeAdded.Count.Indel[j])
		}
	}
}

// findMatchingIndel searches for queryIndel in subjectSlice returning the matching indel in subjectSlice if present, else returning nil if not present
func findMatchingIndel(queryIndel *Indel, subjectSlice []Indel) *Indel {
	var i int
	for i = 0; i < len(subjectSlice); i++ {
		if dna.CompareSeqsIgnoreCase(queryIndel.Alt, subjectSlice[i].Alt) == 0 &&
			dna.CompareSeqsIgnoreCase(queryIndel.Ref, subjectSlice[i].Ref) == 0 {
			return &subjectSlice[i]
		}
	}
	return nil
}

// alleleToVcf converts a *Allele to a vcf record
func alleleToVcf(allele *Allele, p float64, altBase dna.Base, indelSlicePos int, warnings []string) *vcf.Vcf {
	var answer *vcf.Vcf

	filename := path.Base(allele.Sample)
	ext := filepath.Ext(filename)
	name := filename[0 : len(filename)-len(ext)]

	answer = &vcf.Vcf{
		Chr:    allele.Location.Chr,
		Pos:    allele.Location.Pos,
		Id:     ".",
		Qual:   0,
		Filter: delimitStringSlice(warnings, ":"),
		Info:   fmt.Sprintf("Sample=%s:p=%.2e", name, p),
		Format: fmt.Sprintf("DP:AC")}

	switch altBase {
	case dna.A:
		answer.Ref = dna.BaseToString(allele.Count.Ref)
		answer.Alt = "A"
		answer.Notes = fmt.Sprintf("%d:%d", allele.Count.Counts, allele.Count.BaseAF+allele.Count.BaseAR)
	case dna.C:
		answer.Ref = dna.BaseToString(allele.Count.Ref)
		answer.Alt = "C"
		answer.Notes = fmt.Sprintf("%d:%d", allele.Count.Counts, allele.Count.BaseCF+allele.Count.BaseCR)
	case dna.G:
		answer.Ref = dna.BaseToString(allele.Count.Ref)
		answer.Alt = "G"
		answer.Notes = fmt.Sprintf("%d:%d", allele.Count.Counts, allele.Count.BaseGF+allele.Count.BaseGR)
	case dna.T:
		answer.Ref = dna.BaseToString(allele.Count.Ref)
		answer.Alt = "T"
		answer.Notes = fmt.Sprintf("%d:%d", allele.Count.Counts, allele.Count.BaseTF+allele.Count.BaseTR)
	case dna.Gap:
		answer.Ref = dna.BasesToString(allele.Count.Indel[indelSlicePos].Ref)
		answer.Alt = dna.BasesToString(allele.Count.Indel[indelSlicePos].Alt)
		answer.Notes = fmt.Sprintf("%d:%d", allele.Count.Counts, allele.Count.Indel[indelSlicePos].CountF+allele.Count.Indel[indelSlicePos].CountR)
	}
	return answer
}

// alleleFishersExact test performs fishers exact test with the given parameters and returns p
func alleleFishersExact(a int, b int, c int, d int, afThreshold float64) float64 {
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

	// If the allele frequency is less than the threshold then p is noted as 1 so as to be excluded
	case float64(c)/float64(c+a) < afThreshold:
		p = 1

	// If no exclusion conditions are met, then calculate p value
	default:
		p = numbers.FisherExact(a, b, c, d, true)
	}
	return p
}

// getFishersInput determines the variables to use as input to alleleFishersExact
func getFishersInput(experimental *Allele, bkgd *Allele, normalPresent bool, altBase dna.Base, indelSlicePos int) (int, int, int, int) {
	// Begin gathering parameters for Fishers Exact Test done in the numbers package
	// test is for the matrix:
	// [a b]
	// [c d]
	// a = Samples Ref Allele Count
	// b = Background Ref Allele Count - Samples Ref Allele Count
	// c = Samples Alt Allele Count
	// d = Background Alt Allele Count - Samples Alt Allele Count
	var a, b, c, d int

	switch bkgd.Count.Ref {
	case dna.A:
		a = int(experimental.Count.BaseAF + experimental.Count.BaseAR)
		b = int(bkgd.Count.BaseAF + bkgd.Count.BaseAR)
	case dna.C:
		a = int(experimental.Count.BaseCF + experimental.Count.BaseCR)
		b = int(bkgd.Count.BaseCF + bkgd.Count.BaseCR)
	case dna.G:
		a = int(experimental.Count.BaseGF + experimental.Count.BaseGR)
		b = int(bkgd.Count.BaseGF + bkgd.Count.BaseGR)
	case dna.T:
		a = int(experimental.Count.BaseTF + experimental.Count.BaseTR)
		b = int(bkgd.Count.BaseTF + bkgd.Count.BaseTR)
	}

	switch altBase {
	case dna.A:
		c = int(experimental.Count.BaseAF + experimental.Count.BaseAR)
		d = int(bkgd.Count.BaseAF + bkgd.Count.BaseAR)
	case dna.C:
		c = int(experimental.Count.BaseCF + experimental.Count.BaseCR)
		d = int(bkgd.Count.BaseCF + bkgd.Count.BaseCR)
	case dna.G:
		c = int(experimental.Count.BaseGF + experimental.Count.BaseGR)
		d = int(bkgd.Count.BaseGF + bkgd.Count.BaseGR)
	case dna.T:
		c = int(experimental.Count.BaseTF + experimental.Count.BaseTR)
		d = int(bkgd.Count.BaseTF + bkgd.Count.BaseTR)
	case dna.Gap:
		c = int(experimental.Count.Indel[indelSlicePos].CountF + experimental.Count.Indel[indelSlicePos].CountR)
		bkgdIndel := findMatchingIndel(&experimental.Count.Indel[indelSlicePos], bkgd.Count.Indel)
		d = int(bkgdIndel.CountF + bkgdIndel.CountR)
	}

	if !normalPresent {
		b -= a
		d -= c
	}

	return a, b, c, d
}

// delimitStringSlice is a helper function that concatenates a []string with a given delimiter
func delimitStringSlice(strings []string, delimiter string) string {
	var answer string
	var i int = 0
	if len(strings) > 1 {
		for i = 1; i < len(strings); i++ {
			answer = answer + strings[i-1] + delimiter
		}
		answer = answer + strings[i-1]
	} else if len(strings) == 1 {
		answer = answer + strings[0]
	}

	return answer
}

// passesStrandBias returns false if the distribution of forward and reverse reads is more extreme than 95% fwd/rev
func passesStrandBias(alpha int32, beta int32) bool {
	val := float64(alpha) / float64(alpha+beta)
	if val < 0.95 && val > 0.05 {
		return true
	} else {
		return false
	}
}
