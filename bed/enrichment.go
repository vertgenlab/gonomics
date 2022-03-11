package bed

import (
	"github.com/vertgenlab/gonomics/numbers"
	"github.com/vertgenlab/gonomics/numbers/logspace"
	"log"
	"math"
	"strings"
)

// subtractFromBedCoord is a helper function of overlapProbability that subtracts the values subStart and subEnd from entries in a slice of Bed b.
// bClone is a clone of b used to save on memory allocation.
func subtractFromBedCoord(b []Bed, subStart int, subEnd int, bClone []Bed) {
	var prevEnd int = 0
	var prevChrom string = ""

	for i := range b {
		if strings.Compare(prevChrom, "") == 0 || strings.Compare(prevChrom, bClone[i].Chrom) != 0 {
			prevChrom = bClone[i].Chrom
			prevEnd = 0
		}
		bClone[i].Chrom = b[i].Chrom
		bClone[i].ChromStart = numbers.Max(prevEnd, b[i].ChromStart-subStart)
		bClone[i].ChromEnd = numbers.Max(b[i].ChromStart, b[i].ChromEnd-subEnd)
		prevEnd = bClone[i].ChromEnd
	}
}

// overlapProbability calculates the probability that an element of len 'length' overlaps a set of genomic 'elements' in a genome represented by 'noGapRegions'.
// tempElements and tempNoGap represent cloned slices of elements and noGapRegions which are used for memory optimization.
func overlapProbability(elements []Bed, tempElements []Bed, length int, noGapRegions []Bed, tempNoGap []Bed) float64 {
	subtractFromBedCoord(elements, length-1, 0, tempElements)
	subtractFromBedCoord(noGapRegions, 0, length-1, tempNoGap)
	return float64(OverlapLengthSum(tempElements, tempNoGap)) / float64(TotalSize(tempNoGap))
}

// ElementOverlapProbabilities returns a slice of float64 representing the probabilities that an element in elements2 overlaps an element in elements1.
func ElementOverlapProbabilities(elements1 []Bed, elements2 []Bed, noGapRegions []Bed) []float64 {
	var answer []float64 = make([]float64, len(elements2))
	var tempElements1 []Bed = make([]Bed, len(elements1))
	copy(tempElements1, elements1)
	var tempNoGap []Bed = make([]Bed, len(noGapRegions))
	copy(tempNoGap, noGapRegions)
	var tempElements2 []Bed = make([]Bed, len(elements2))
	copy(tempElements2, elements2)
	SortBySize(tempElements2)
	var currLen, prevLen int = 0, 0

	for i := range tempElements2 {
		currLen = tempElements2[i].ChromEnd - tempElements2[i].ChromStart
		if currLen == prevLen {
			answer[i] = answer[i-1]
		} else {
			answer[i] = overlapProbability(elements1, tempElements1, currLen, noGapRegions, tempNoGap)
			prevLen = currLen
		}
	}

	return answer
}

// EnrichmentPValueApproximation performs a calculation of enrichment based on a set of overlap probabilities and the observed overlap count.
// Uses a normal approximation of the resulting binomial distribution. Very fast and should converge on the correct answer. This is the preferred method.
// Returns a slice of three values. The first represents a debug check (hardcoded to one), the second is the expected number of overlaps, and the third is the pValue of the observed overlap.
func EnrichmentPValueApproximation(elementOverlapProbs []float64, overlapCount int) []float64 {
	var answer []float64 = make([]float64, 3)
	var mu, sigma float64 = 0, 0

	for i := range elementOverlapProbs {
		mu += elementOverlapProbs[i]
		sigma = sigma + elementOverlapProbs[i]*(1-elementOverlapProbs[i])
	}

	sigma = math.Sqrt(sigma) //sigma will represent the standard deviation for our normal approximation. the above sum gives us the variance.

	answer[0] = 1.0
	answer[1] = mu //mu represents the expected value

	//calculate pValue approximation
	pValue := numbers.NormalDist(float64(overlapCount), mu, sigma)
	for s := overlapCount + 1; s <= len(elementOverlapProbs); s++ {
		pValue += numbers.NormalDist(float64(s), mu, sigma)
	}
	answer[2] = pValue

	return answer
}

// EnrichmentPValueExact performs an exact calculation fo enrichment bgaased on a set of overlap probabilities and the observed overlap count.
// The exact method is non-polynomial, and is thus not recommended for large datasets.
// Returns a slice of three values. The first is the debug check, the second is the expected number of overlaps, and the third is the pValue of the observed number of overlaps.
func EnrichmentPValueExact(elementOverlapProbs []float64, overlapCount int) []float64 {
	var numTrials = len(elementOverlapProbs)
	var answer []float64 = make([]float64, 3)
	var prevCol []float64 = make([]float64, numTrials+1)
	var currCol []float64 = make([]float64, numTrials+1)

	//for one trial (a bed input with a single element)
	prevCol[0] = math.Log(1 - elementOverlapProbs[0])
	currCol[0] = prevCol[0]
	prevCol[1] = math.Log(elementOverlapProbs[0])
	currCol[1] = math.Log(elementOverlapProbs[0])

	var s int
	var expected, pValue float64

	for t := 1; t < numTrials; t++ {
		prevCol, currCol = currCol, prevCol
		currCol[0] = prevCol[0] + math.Log(1-elementOverlapProbs[t])
		for s = 1; s <= t; s++ {
			currCol[s] = logspace.Add(prevCol[s]+math.Log(1-elementOverlapProbs[t]), prevCol[s-1]+math.Log(elementOverlapProbs[t]))
		}
		currCol[t+1] = prevCol[t] + math.Log(elementOverlapProbs[t])
	}

	//Now we check whether our probabilities add to 1.
	check := currCol[0]
	for s = 1; s <= numTrials; s++ {
		check = logspace.Add(check, currCol[s])
		if s == 1 {
			expected = currCol[s]
		} else {
			expected = logspace.Add(expected, currCol[s]+math.Log(float64(s)))
		}
	}

	pValue = currCol[overlapCount]
	for s = overlapCount + 1; s <= numTrials; s++ {
		pValue = logspace.Add(pValue, currCol[s])
	}

	answer[0] = math.Exp(check)
	answer[2] = math.Exp(pValue)
	answer[1] = math.Exp(expected)
	return answer
}

// EnrichmentPValueUpperBound, together with EnrichmentPValueLowerBound, provide a range of possible values for the pValue of overlap.
// Returns a slice of three values. The first is the debug check, the second is the expected number of overlaps, and the third is the pValue of the observed number of overlaps.
func EnrichmentPValueUpperBound(elements1 []Bed, elements2 []Bed, noGapRegions []Bed, overlapCount int, verbose int) []float64 {
	var numTrials int = len(elements2)
	var answer []float64 = make([]float64, 3)
	var tempElements1 []Bed = make([]Bed, len(elements1))
	copy(tempElements1, elements1)
	var tempNoGap []Bed = make([]Bed, len(noGapRegions))
	copy(tempNoGap, noGapRegions)
	minElements2 := findLargestBedLength(elements2)
	var curr float64

	if verbose > 0 {
		log.Println("Calculating overlapProbability.")
	}
	prob := overlapProbability(elements1, tempElements1, minElements2, noGapRegions, tempNoGap)

	if verbose > 0 {
		log.Println("Calculating the pValue.")
	}

	pValue, underflow := numbers.BinomialDist(numTrials, overlapCount, prob)
	for s := overlapCount + 1; s <= numTrials && !(underflow && float64(s) > float64(numTrials)*prob); s++ {
		curr, underflow = numbers.BinomialDist(numTrials, s, prob)
		pValue += curr
	}

	answer[0] = 1 //hardcoded for now, we don't do the check with this method.
	answer[2] = pValue
	answer[1] = prob * float64(numTrials)
	return answer
}

// EnrichmentPValueLowerBound, together with EnrichmentPValueUpperBound, provide a range of possible values for the pValue of overlap.
// Returns a slice of three values. The first is the debug check, the second is the expected number of overlaps, and the third is the pValue of the observed number of overlaps.
func EnrichmentPValueLowerBound(elements1 []Bed, elements2 []Bed, noGapRegions []Bed, overlapCount int, verbose int) []float64 {
	var numTrials int = len(elements2)
	var answer []float64 = make([]float64, 3)
	var tempElements1 []Bed = make([]Bed, len(elements1))
	var tempNoGap []Bed = make([]Bed, len(noGapRegions))
	copy(tempElements1, elements1)
	copy(tempNoGap, noGapRegions)
	minElements2 := findShortestBedLength(elements2)
	if verbose > 0 {
		log.Println("Calculating overlapProbability.")
	}
	prob := overlapProbability(elements1, tempElements1, minElements2, noGapRegions, tempNoGap)

	if verbose > 0 {
		log.Println("Calculating the pValue.")
	}
	var curr float64
	pValue, underflow := numbers.BinomialDist(numTrials, overlapCount, prob)
	if underflow {
		answer[2] = 0.0
		if verbose > 0 {
			log.Println("Underflow detected in the binomial distribution at overlapCount. p value is too small to detect.")
		}
	} else {
		for s := overlapCount + 1; s <= numTrials; s++ {
			curr, underflow = numbers.BinomialDist(numTrials, s, prob)
			if underflow {
				break
			}
			pValue += curr
		}
	}

	answer[0] = 1 //hardcoded for now, we don't do the check with this method.
	answer[2] = pValue
	answer[1] = prob * float64(numTrials)
	return answer
}

// helperfunction that returns an int corresponding to the longest distance between chromstart and chromEnd from an input set of bed entries.
func findLargestBedLength(b []Bed) int {
	var maxLength int = b[0].ChromEnd - b[0].ChromStart

	for i := 1; i < len(b); i++ {
		if b[i].ChromEnd-b[i].ChromStart > maxLength {
			maxLength = b[i].ChromEnd - b[i].ChromStart
		}
	}
	return maxLength
}

// helperfunction that returns an int corresponding to the shortest distance between chromstart and chromEnd from an input set of bed entries.
func findShortestBedLength(b []Bed) int {
	var minLength int = b[0].ChromEnd - b[0].ChromStart
	for i := 1; i < len(b); i++ {
		if b[i].ChromEnd-b[i].ChromStart < minLength {
			minLength = b[i].ChromEnd - b[i].ChromStart
		}
	}
	return minLength
}
