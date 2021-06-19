package bed

import (
	"github.com/vertgenlab/gonomics/numbers"
	"math"
	"strings"
)

//subtractFromBedCoord is a helper function of overlapProbability that subtracts the values subStart and subEnd from entries in a slice of Bed b.
//bClone is a clone of b used to save on memory allocation.
func subtractFromBedCoord(b []*Bed, subStart int, subEnd int, bClone []*Bed) {
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

//overlapProbability calculates the probability that an element of len 'length' overlaps a set of genomic 'elements' in a genome represented by 'noGapRegions'.
//tempElements and tempNoGap represent cloned slices of elements and noGapRegions which are used for memory optimization.
func overlapProbability(elements []*Bed, tempElements []*Bed, length int, noGapRegions []*Bed, tempNoGap []*Bed) float64 {
	subtractFromBedCoord(elements, length-1, 0, tempElements)
	subtractFromBedCoord(noGapRegions, 0, length-1, tempNoGap)
	return float64(OverlapLengthSum(tempElements, tempNoGap)) / float64(TotalSize(tempNoGap))
}

func ElementOverlapProbabilities(elements1 []*Bed, elements2 []*Bed, noGapRegions []*Bed) []float64 {
	var tempElements1 []*Bed = make([]*Bed, len(elements1))
	var tempElements2 []*Bed = make([]*Bed, len(elements2))
	var tempNoGap []*Bed = make([]*Bed, len(noGapRegions))
	var answer []float64 = make([]float64, len(elements2))
	copy(tempElements1, elements1)
	copy(tempNoGap, noGapRegions)
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

func EnrichmentPValue(elementOverlapProbs []float64, numTrials int, overlapCount int) []float64 {
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
			currCol[s] = numbers.AddLog(prevCol[s]+math.Log(1-elementOverlapProbs[t]), prevCol[s-1]+math.Log(elementOverlapProbs[t]))
		}
		currCol[t+1] = prevCol[t] + math.Log(elementOverlapProbs[t])
	}

	//Now we check whether our probabilities add to 1.
	check := currCol[0]
	for s = 1; s <= numTrials; s++ {
		check = numbers.AddLog(check, currCol[s])
		if s == 1 {
			expected = currCol[s]
		} else {
			expected = numbers.AddLog(expected, currCol[s]+math.Log(float64(s)))
		}
	}

	pValue = currCol[overlapCount]
	for s = overlapCount + 1; s <= numTrials; s++ {
		pValue = numbers.AddLog(pValue, currCol[s])
	}

	answer[0] = math.Exp(check)
	answer[2] = math.Exp(pValue)
	answer[1] = math.Exp(expected)
	return answer
}
