// Package simulate contains functions for simulation of genomic data with ils
// variants, reads, or regions.
package simulate

import (
	"log"
	"math/rand"

	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/expandedTree"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/genePred"
	"github.com/vertgenlab/gonomics/numbers"
)

func sumArray(numbers []int) int {
    result := 0
    for i := 0; i < len(numbers); i++ {
        result += numbers[i]
    }
    return result
}

// takes as input
// a list of topologies (v0, ...)
// desired length
// desired proportion of each topology
// need to calculate the total number of bases in each topology (according to the provided proportions)
// and then for each topology's proportional length randomly generate a list of numbers that sum to that length
// and also that the total length of that list of numbers (ie the number of numbers in that list), over all topologies is equal to numSegs
func SimulateIls(roots []*expandedTree.ETree, proportions []float64, numSegs int, totalLength int, seed int, GC float64) {
	rand.Seed(seed)
	/// need to a) generate numSegs number of segments that in total sum to totalLength 
	
	if len(roots) != len(proportions) {
		log.Fatal("Must provide proportions same length as roots, can be zero")
	}

	numTopos := len(roots)

	if sumArray(proportions) != 1 {
		log.Fatal("Must provide proportions that sum to 1")
	}

	if totalLength < 1 {
		log.Fatal("Cannot generate sequence with non-positive length")
	} 
	// put default value in cli
	
	// else if totalLength == 0 {
	// 	// randomly generate total length > 1000000
	// 	totalLength = rand.Intn(5000000) + 1000000
	// }

	if numSegs < 0 {
		log.Fatal("Cannot request less than 1 segments")
	}
	// for practical applications
	// would be good not to have random unreported-to-user behaviour
	// put this in cli
	
	// else if numSegs == 0 {
	// 	// randomly generate number of segments
	// 	// tries to assume that segments are length 3bps each
	// 	numSegs = rand.Intn(totalLength - totalLength/3 + 1) + totalLength/3
	// }

	/// calculate the total number of bases in each topology
	topoLength := make([]int, numTopos)
	remainders := make([]float64, numTopos)
	sum := 0

	/// calculate how many bases will be in each topo, needs to be integer
	// e.g. proportions for testcase will be length 4 for V0/1/2/3 each
	for i, fraction := range proportions {
		fractionLength := fraction * totalLength
		topoLength[i] = int(fractionLength) // floor
		remainders[i] = fractionLength - float64(topoLength[i])
		sum += topoLength[i]
	}

	/// distribute the remaining if doesn't sum to int, based on the largest remainder
	// at most missing is 3
	missing := totalLength - sum
	for missing > 0 {
		maxIdx := 0
		for i := 1; i < numTopos; i++ {
			if remainders[i] > remainders[maxIdx] {
				maxIdx = i
			}
		}
		topoLength[maxIdx]++ // increment the topology's length by 1
		remainders[maxIdx] = 0 // avoid picking again
		missing--
		// keep on picking the larget remainder until missing==0
	}


	/// randomly calculate the number of segments in each topology, must sum to numSegs over all topologies
	/// at most will allow numSegment per topo == corresponding length in topoLength
	topoLengthdivider = 1
	tempLength := totalLength
	// if topoLength >= 10, topoLengthDivisor is 10
	for tempLength > 99 {
		tempLength /= 10
		topoLengthdivider *= 10
	}

	numSegsPerTopo := make([]int, numTopos)
	for idx := numSegsPerTopo {
		numSegsPerTopo[idx] = topoLength[idx]/topoLengthDivider ????? need to randomise this
	}

	remainingSegs = numSegs - topoLengthDivider*4




	centerDistrib := make([]float64, numTopos)
	for i := range topoLength {
		centerDistrib[i] = float64(numSegs) * float64(topoLength[i]) / float64(totalCap)
	}

	// Distribute the remaining segments
	for remaining > 0 {
		totalWeight := 0.0
		weights := make([]float64, numTopos)

		for i := 0; i < numTopos; i++ {
			if numSegsPerTopo[i] >= topoLength[i] {
				weights[i] = 0
				continue
			}

			// How far below target are we?
			deficit := centerDistrib[i] - float64(numSegsPerTopo[i])

			// Bias toward under-target topologies, but still allow randomness.
			// Increase/decrease the constants to make it tighter/looser.
			w := 0.25
			if deficit > 0 {
				w += deficit
			}

			weights[i] = w
			totalWeight += w
		}

		if totalWeight == 0 {
			panic("no available topology has spare capacity")
		}

		r := rand.Float64() * totalWeight
		for i, w := range weights {
			r -= w
			if r <= 0 {
				numSegsPerTopo[i]++
				remaining--
				break
			}
		}
	}

	inputSeqs := make([]fasta.Fasta, numTopos)
	for _, lenSeq in topoLength {
		fasta.Fasta{Name: fmt.Sprintf("Sequence_%v", i), Seq: simulate.RandIntergenicSeq(GC, lenSeq)}
	}




}