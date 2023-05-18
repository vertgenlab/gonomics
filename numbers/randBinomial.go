package numbers

import (
	"math/rand"
)

// BinomialAlias contains alias table information used to
// generate binomial distributed random variates for pre-specified
// parameters. Setting up the BinomialAlias struct runs in O(n) time.
// However, once the alias is set up for a fixed cost, random variates
// can be generated from the distribution in O(1) time.
// More information on the alias method can be found here:
// https://en.wikipedia.org/wiki/Alias_method
type BinomialAlias struct {
	Probability []float64
	Alias       []int
}

// RandBinomial generates binomial distributed variates from
// a pre-generated BinomialAlias struct, which can be made for
// a specified binomial distribution with 'MakeBinomialAlias'.
func RandBinomial(alias BinomialAlias) int {
	index := RandIntInRange(0, len(alias.Probability))
	if rand.Float64() < alias.Probability[index] {
		return index
	} else {
		return alias.Alias[index]
	}
}

// MakeBinomialAlias generates a BinomialAlias struct for a specified binomial
// distribution of n trials with success probability p.
// Note if the probability of a binomial outcome is less than math.MinFloat64, it will never be generated.
// This implementation therefore draws variates from an approximation of the binomial distribution
// that is truncated to zero when its tail crosses below the float underflow threshold.
func MakeBinomialAlias(n int, p float64) BinomialAlias {
	var underflow bool
	var currIndex int
	var currOver, currUnder int
	var emptyRoomInUnderFull float64

	answer := BinomialAlias{Probability: make([]float64, n+1), Alias: make([]int, n+1)} // n + 1 possible outcomes for n trials
	oneOverNPlusOne := 1.0 / float64(n+1)                                               // oneOverNPlusOne is the amount of probability in a full bucket

	for currIndex = range answer.Probability {
		answer.Probability[currIndex], underflow = BinomialDist(n, currIndex, p)
		if underflow {
			answer.Probability[currIndex] = 0
		}
	}

	//two stacks are made. One for probabilities greater than average (oneOverN) and the other for probabilities less than or equal to average
	var underFull = make([]int, 0)
	var overFull = make([]int, 0)
	for currIndex = range answer.Probability {
		if answer.Probability[currIndex] > oneOverNPlusOne {
			overFull = append(overFull, currIndex)
		} else {
			underFull = append(underFull, currIndex)
		}
	}

	for len(overFull) > 0 && len(underFull) > 0 {
		currUnder = underFull[len(underFull)-1] // last entry in underFull
		currOver = overFull[len(overFull)-1]    // last entry in overFull
		emptyRoomInUnderFull = oneOverNPlusOne - answer.Probability[currUnder]

		//pour last overFull into last underFull, filling the first underFull
		answer.Alias[currUnder] = currOver       // set alias of answer (in other words, fill the bucket)
		underFull = underFull[:len(underFull)-1] // remove last slice entry
		answer.Probability[currOver] -= emptyRoomInUnderFull

		//if currOver is now underFull, move to that slice
		if answer.Probability[currOver] < oneOverNPlusOne {
			underFull = append(underFull, currOver)
			overFull = overFull[:len(overFull)-1]
		}
	}

	// should either overFull or underFull still have entries (like when
	// the starting probability is 1/(N+1)), we set the probability to
	// 1/(n+1) and remove from the list.
	for len(overFull) > 0 {
		currOver = overFull[len(overFull)-1]
		answer.Probability[currOver] = oneOverNPlusOne
		overFull = overFull[:len(overFull)-1]
	}
	for len(underFull) > 0 {
		currUnder = underFull[len(underFull)-1]
		answer.Probability[currUnder] = oneOverNPlusOne
		underFull = underFull[:len(underFull)-1]
	}

	//We now normalize the probabilities so each bucket ranges from 0 to 1
	// instead of 0 to 1/(N+1)
	for currIndex = range answer.Probability {
		answer.Probability[currIndex] = answer.Probability[currIndex] * float64(n+1.0)
	}

	return answer
}
