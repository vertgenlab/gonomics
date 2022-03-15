package lift

import (
	"github.com/vertgenlab/gonomics/numbers"
	"github.com/vertgenlab/gonomics/numbers/logspace"
	"log"
	"math"
	"strings"
)

// ElementOverlapProbabilities returns a slice of float64 representing the probabilities that an element in elements2 overlaps an element in elements1.
func ElementOverlapProbabilities(elements1 []Lift, elements2 []Lift, noGapRegions []Lift) []float64 {
	var answer []float64 = make([]float64, len(elements2))
	var tempElements1 []Lift = make([]Lift, len(elements1))
	copy(tempElements1, elements1)
	var tempNoGap []Lift = make([]Lift, len(noGapRegions))
	copy(tempNoGap, noGapRegions)
	var tempElements2 []Lift = make([]Lift, len(elements2))
	copy(tempElements2, elements2)
	sortBySize(tempElements2)
	var currLen, prevLen int = 0, 0

	for i := range tempElements2 {
		currLen = tempElements2[i].GetChromEnd() - tempElements2[i].GetChromStart()
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
// Returns a slice of four values. The first represents a debug check (hardcoded to one), the second is the expected number of overlaps, and the third and fourth represent the pValues for enrichment and depletion, respectively.
func EnrichmentPValueApproximation(elementOverlapProbs []float64, overlapCount int) []float64 {
	var answer []float64 = make([]float64, 4)
	var mu, sigma float64 = 0, 0
	var enrichPValue, depletePValue float64
	var s int

	for i := range elementOverlapProbs {
		mu += elementOverlapProbs[i]
		sigma = sigma + elementOverlapProbs[i]*(1-elementOverlapProbs[i])
	}

	sigma = math.Sqrt(sigma) //sigma will represent the standard deviation for our normal approximation. the above sum gives us the variance.

	answer[0] = 1.0
	answer[1] = mu //mu represents the expected value

	//calculate pValue approximation
	enrichPValue = numbers.NormalDist(float64(overlapCount), mu, sigma)
	for s = overlapCount + 1; s <= len(elementOverlapProbs); s++ {
		enrichPValue += numbers.NormalDist(float64(s), mu, sigma)
	}
	answer[2] = enrichPValue

	depletePValue = numbers.NormalDist(float64(overlapCount), mu, sigma)
	for s = overlapCount - 1; s >= 0; s-- {
		depletePValue += numbers.NormalDist(float64(s), mu, sigma)
	}
	answer[3] = depletePValue

	return answer
}

// EnrichmentPValueExact performs an exact calculation fo enrichment bgaased on a set of overlap probabilities and the observed overlap count.
// The exact method is non-polynomial, and is thus not recommended for large datasets.
// Returns a slice of four values. The first is the debug check, the second is the expected number of overlaps, and the third and fourth represent the pValues for enrichment and depletion, respectively.
func EnrichmentPValueExact(elementOverlapProbs []float64, overlapCount int) []float64 {
	var numTrials = len(elementOverlapProbs)
	var answer []float64 = make([]float64, 4)
	var prevCol []float64 = make([]float64, numTrials+1)
	var currCol []float64 = make([]float64, numTrials+1)

	//for one trial (a bed input with a single element)
	prevCol[0] = math.Log(1 - elementOverlapProbs[0])
	currCol[0] = prevCol[0]
	prevCol[1] = math.Log(elementOverlapProbs[0])
	currCol[1] = math.Log(elementOverlapProbs[0])

	var s int
	var expected, enrichPValue, depletePValue float64

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

	enrichPValue = currCol[overlapCount]
	for s = overlapCount + 1; s <= numTrials; s++ {
		enrichPValue = logspace.Add(enrichPValue, currCol[s])
	}

	depletePValue = currCol[overlapCount]
	for s = overlapCount - 1; s >= 0; s-- {
		depletePValue = logspace.Add(depletePValue, currCol[s])
	}

	answer[0] = math.Exp(check)
	answer[1] = math.Exp(expected)
	answer[2] = math.Exp(enrichPValue)
	answer[3] = math.Exp(depletePValue)

	return answer
}

// EnrichmentPValueUpperBound, together with EnrichmentPValueLowerBound, provide a range of possible values for the pValue of overlap.
// Returns a slice of four values. The first is the debug check, the second is the expected number of overlaps, and the third and fourth represent the pValues for enrichment and depletion, respectively.
func EnrichmentPValueUpperBound(elements1 []Lift, elements2 []Lift, noGapRegions []Lift, overlapCount int, verbose int) []float64 {
	var numTrials int = len(elements2)
	var answer []float64 = make([]float64, 4)
	var tempElements1 []Lift = make([]Lift, len(elements1))
	var tempNoGap []Lift = make([]Lift, len(noGapRegions))
	var curr, enrichPValue, depletePValue float64
	var s int

	copy(tempElements1, elements1)
	copy(tempNoGap, noGapRegions)

	minElements2 := findLargestLength(elements2)

	if verbose > 0 {
		log.Println("Calculating overlapProbability.")
	}
	prob := overlapProbability(elements1, tempElements1, minElements2, noGapRegions, tempNoGap)

	if verbose > 0 {
		log.Println("Calculating the pValue.")
	}

	enrichPValue, _ = numbers.BinomialDist(numTrials, overlapCount, prob)
	for s = overlapCount + 1; s <= numTrials; s++ {
		curr, _ = numbers.BinomialDist(numTrials, s, prob)
		enrichPValue += curr
	}

	depletePValue, _ = numbers.BinomialDist(numTrials, overlapCount, prob)
	for s = overlapCount - 1; s >= 0; s-- {
		curr, _ = numbers.BinomialDist(numTrials, s, prob)
		depletePValue += curr
	}

	answer[0] = 1 //hardcoded for now, we don't do the check with this method.
	answer[1] = prob * float64(numTrials)
	answer[2] = enrichPValue
	answer[3] = depletePValue
	return answer
}

// EnrichmentPValueLowerBound, together with EnrichmentPValueUpperBound, provide a range of possible values for the pValue of overlap.
// Returns a slice of four values. The first is the debug check, the second is the expected number of overlaps, and the third and fourth represent the pValues for enrichment and depletion, respectively.
func EnrichmentPValueLowerBound(elements1 []Lift, elements2 []Lift, noGapRegions []Lift, overlapCount int, verbose int) []float64 {
	var numTrials int = len(elements2)
	var answer []float64 = make([]float64, 4)
	var tempElements1 []Lift = make([]Lift, len(elements1))
	var tempNoGap []Lift = make([]Lift, len(noGapRegions))
	var enrichPValue, depletePValue float64
	var s int

	copy(tempElements1, elements1)
	copy(tempNoGap, noGapRegions)

	minElements2 := findShortestLength(elements2)

	if verbose > 0 {
		log.Println("Calculating overlapProbability.")
	}
	prob := overlapProbability(elements1, tempElements1, minElements2, noGapRegions, tempNoGap)

	if verbose > 0 {
		log.Println("Calculating the pValue.")
	}
	var curr float64
	enrichPValue, _ = numbers.BinomialDist(numTrials, overlapCount, prob)
	for s = overlapCount + 1; s <= numTrials; s++ {
		curr, _ = numbers.BinomialDist(numTrials, s, prob)
		enrichPValue += curr
	}

	depletePValue, _ = numbers.BinomialDist(numTrials, overlapCount, prob)
	for s = overlapCount - 1; s >= 0; s-- {
		curr, _ = numbers.BinomialDist(numTrials, s, prob)
		depletePValue += curr
	}

	answer[0] = 1 //hardcoded for now, we don't do the check with this method.
	answer[1] = prob * float64(numTrials)
	answer[2] = enrichPValue
	answer[3] = depletePValue

	return answer
}

func subtractFromCoord(regions []Lift, subStart int, subEnd int, regionClone []Lift) {
	var prevEnd int = 0
	var prevChrom string = ""

	for i := range regions {
		if strings.Compare(prevChrom, "") == 0 || strings.Compare(prevChrom, regionClone[i].GetChrom()) != 0 {
			prevChrom = regionClone[i].GetChrom()
			prevEnd = 0
		}
		regionClone[i] = regionClone[i].UpdateCoord(regions[i].GetChrom(), numbers.Max(prevEnd, regions[i].GetChromStart()-subStart), numbers.Max(regions[i].GetChromStart(), regions[i].GetChromEnd()-subEnd)).(Lift)
		prevEnd = regionClone[i].GetChromEnd()
	}
}
