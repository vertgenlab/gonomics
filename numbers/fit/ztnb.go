package fit

import (
	"fmt"
	"github.com/vertgenlab/gonomics/numbers"
	"github.com/vertgenlab/gonomics/numbers/logspace"
	"log"
	"math"
	"math/rand"
)

// direction is a type used to define the heading in coordinate ascent.
type direction byte

const (
	neutral   direction = 0
	north     direction = 1
	northeast direction = 2
	east      direction = 3
	southeast direction = 4
	south     direction = 5
	southwest direction = 6
	west      direction = 7
	northwest direction = 8
)

// plotLossSurfaceZTNB plots a matrix of log-likelihood values for a ZTNB fit to input histogram format data from the range
// rMax to rMin with resolution rStep and from pMin to pMax with resolution pStep. Useful for debugging.
func plotLossSurfaceZTNB(data []int, rMin float64, rMax float64, rStep float64, pMin float64, pMax float64, pStep float64) ([][]float64, float64, float64, float64) {
	var answer [][]float64 = make([][]float64, int((pMax-pMin)/pStep))
	lowestLoss := zeroTruncatedNegativeBinomialLogLikelihood(data, rMin, pMin)
	lowestLossR := rMin
	lowestLossP := pMin
	var currLoss float64
	for i := range answer {
		answer[i] = make([]float64, int((rMax-rMin)/rStep))
		for j := range answer[i] {
			currLoss = zeroTruncatedNegativeBinomialLogLikelihood(data, rMin+rStep*float64(i), pMin+pStep*float64(j))
			answer[i][j] = currLoss
			if currLoss > lowestLoss {
				lowestLoss = currLoss
				lowestLossR = rMin + rStep*float64(i)
				lowestLossP = pMin + pStep*float64(j)
			}
		}
	}

	return answer, lowestLoss, lowestLossR, lowestLossP
}

// zeroTruncatedNegativeBinomialLogLikelihood returns the natural logarithm of the likelihood of input histogram
// count data given a zero-truncated negative binomial (ZTNB) distribution parameterized by R and P.
func zeroTruncatedNegativeBinomialLogLikelihood(data []int, R float64, P float64) float64 {
	likelihood := 0.0
	var density float64
	for i := 1; i < len(data); i++ {
		density, _ = numbers.NegativeBinomialDist(i, R, P, true)
		if 1-math.Pow(P, R) == 0 {
			fmt.Printf("Found it! P: %v. R: %v.\n", P, R)
		}
		likelihood += float64(data[i]) * logspace.Divide(density, math.Log(1-math.Pow(P, R)))
	}
	return likelihood
}

// checkNorth is a helper function of nextDirection and firstDirection. This function checks the likelihood in the 'north'
// direction, when R is higher than the current parameter set by rStep.
func checkNorth(data []int, R float64, P float64, rStep float64, pStep float64, currDirection direction, currLikelihood float64) (direction, float64) {
	nextLikelihood := zeroTruncatedNegativeBinomialLogLikelihood(data, R+rStep, P)
	if nextLikelihood > currLikelihood {
		return north, nextLikelihood
	}
	return currDirection, currLikelihood
}

// checkNorthEast is a helper function of nextDirection and firstDirection. This function checks the likelihood in the 'northeast'
// direction, when R and P are higher than the current parameter set by rStep and pStep.
func checkNorthEast(data []int, R float64, P float64, rStep float64, pStep float64, currDirection direction, currLikelihood float64) (direction, float64) {
	if P+pStep > 0.999 {
		return currDirection, currLikelihood
	}
	nextLikelihood := zeroTruncatedNegativeBinomialLogLikelihood(data, R+rStep, P+pStep)
	if nextLikelihood > currLikelihood {
		return northeast, nextLikelihood
	}
	return currDirection, currLikelihood
}

// checkEast is a helper function of nextDirection and firstDirection. This function checks the likelihood in the 'east'
// direction, when P is higher than the current parameter set by pStep.
func checkEast(data []int, R float64, P float64, rStep float64, pStep float64, currDirection direction, currLikelihood float64) (direction, float64) {
	if P+pStep > 0.999 {
		return currDirection, currLikelihood
	}
	nextLikelihood := zeroTruncatedNegativeBinomialLogLikelihood(data, R, P+pStep)
	if nextLikelihood > currLikelihood {
		return east, nextLikelihood
	}
	return currDirection, currLikelihood
}

// checkSouthEast is a helper function of nextDirection and firstDirection. This function checks the likelihood in the 'southeast'
// direction, when P is higher by pStep and R is lower by rStep.
func checkSouthEast(data []int, R float64, P float64, rStep float64, pStep float64, currDirection direction, currLikelihood float64) (direction, float64) {
	if P+pStep > 0.999 || R-rStep < 0.001 {
		return currDirection, currLikelihood
	}
	nextLikelihood := zeroTruncatedNegativeBinomialLogLikelihood(data, R-rStep, P+pStep)
	if nextLikelihood > currLikelihood {
		return southeast, nextLikelihood
	}
	return currDirection, currLikelihood
}

// checkSouth is a helper function of nextDirection and firstDirection. This function checks the likelihood in the 'south'
// direction, when R is lower by rStep.
func checkSouth(data []int, R float64, P float64, rStep float64, pStep float64, currDirection direction, currLikelihood float64) (direction, float64) {
	if R-rStep < 0.001 {
		return currDirection, currLikelihood
	}
	nextLikelihood := zeroTruncatedNegativeBinomialLogLikelihood(data, R-rStep, P)
	if nextLikelihood > currLikelihood {
		return south, nextLikelihood
	}
	return currDirection, currLikelihood
}

func checkSouthWest(data []int, R float64, P float64, rStep float64, pStep float64, currDirection direction, currLikelihood float64) (direction, float64) {
	if R-rStep < 0.001 || P-pStep < 0.001 {
		return currDirection, currLikelihood
	}
	nextLikelihood := zeroTruncatedNegativeBinomialLogLikelihood(data, R-rStep, P-pStep)
	if nextLikelihood > currLikelihood {
		return southwest, nextLikelihood
	}
	return currDirection, currLikelihood
}

func checkWest(data []int, R float64, P float64, rStep float64, pStep float64, currDirection direction, currLikelihood float64) (direction, float64) {
	if P-pStep < 0.001 {
		return currDirection, currLikelihood
	}
	nextLikelihood := zeroTruncatedNegativeBinomialLogLikelihood(data, R, P-pStep)
	if nextLikelihood > currLikelihood {
		return west, nextLikelihood
	}
	return currDirection, currLikelihood
}

func checkNorthWest(data []int, R float64, P float64, rStep float64, pStep float64, currDirection direction, currLikelihood float64) (direction, float64) {
	if P-pStep < 0.001 {
		return currDirection, currLikelihood
	}
	nextLikelihood := zeroTruncatedNegativeBinomialLogLikelihood(data, R+rStep, P-pStep)
	if nextLikelihood > currLikelihood {
		return northwest, nextLikelihood
	}
	return currDirection, currLikelihood
}

// nextDirection is a helper function of ZeroTruncatedNegativeBinomial. For a given data, parameter set R and P,
// and current direction currDirection, this function returns the optimal direction of travel for the next iteration
// and the likelihood at that destination.
func nextDirection(data []int, R float64, P float64, rStep float64, pStep float64, currDirection direction, currLikelihood float64) (direction, float64) {
	prevLikelihood := currLikelihood
	switch currDirection {
	case neutral:
		return neutral, currLikelihood
	case north:
		currDirection, currLikelihood = checkWest(data, R, P, rStep, pStep, currDirection, currLikelihood)
		currDirection, currLikelihood = checkNorthWest(data, R, P, rStep, pStep, currDirection, currLikelihood)
		currDirection, currLikelihood = checkNorth(data, R, P, rStep, pStep, currDirection, currLikelihood)
		currDirection, currLikelihood = checkNorthEast(data, R, P, rStep, pStep, currDirection, currLikelihood)
		currDirection, currLikelihood = checkEast(data, R, P, rStep, pStep, currDirection, currLikelihood)
	case northwest:
		currDirection, currLikelihood = checkSouthWest(data, R, P, rStep, pStep, currDirection, currLikelihood)
		currDirection, currLikelihood = checkWest(data, R, P, rStep, pStep, currDirection, currLikelihood)
		currDirection, currLikelihood = checkNorthWest(data, R, P, rStep, pStep, currDirection, currLikelihood)
		currDirection, currLikelihood = checkNorth(data, R, P, rStep, pStep, currDirection, currLikelihood)
		currDirection, currLikelihood = checkNorthEast(data, R, P, rStep, pStep, currDirection, currLikelihood)
	case west:
		currDirection, currLikelihood = checkSouth(data, R, P, rStep, pStep, currDirection, currLikelihood)
		currDirection, currLikelihood = checkSouthWest(data, R, P, rStep, pStep, currDirection, currLikelihood)
		currDirection, currLikelihood = checkWest(data, R, P, rStep, pStep, currDirection, currLikelihood)
		currDirection, currLikelihood = checkNorthWest(data, R, P, rStep, pStep, currDirection, currLikelihood)
		currDirection, currLikelihood = checkNorth(data, R, P, rStep, pStep, currDirection, currLikelihood)
	case southwest:
		currDirection, currLikelihood = checkSouthEast(data, R, P, rStep, pStep, currDirection, currLikelihood)
		currDirection, currLikelihood = checkSouth(data, R, P, rStep, pStep, currDirection, currLikelihood)
		currDirection, currLikelihood = checkSouthWest(data, R, P, rStep, pStep, currDirection, currLikelihood)
		currDirection, currLikelihood = checkWest(data, R, P, rStep, pStep, currDirection, currLikelihood)
		currDirection, currLikelihood = checkNorthWest(data, R, P, rStep, pStep, currDirection, currLikelihood)
	case south:
		currDirection, currLikelihood = checkEast(data, R, P, rStep, pStep, currDirection, currLikelihood)
		currDirection, currLikelihood = checkSouthEast(data, R, P, rStep, pStep, currDirection, currLikelihood)
		currDirection, currLikelihood = checkSouth(data, R, P, rStep, pStep, currDirection, currLikelihood)
		currDirection, currLikelihood = checkSouthWest(data, R, P, rStep, pStep, currDirection, currLikelihood)
		currDirection, currLikelihood = checkWest(data, R, P, rStep, pStep, currDirection, currLikelihood)
	case southeast:
		currDirection, currLikelihood = checkNorthEast(data, R, P, rStep, pStep, currDirection, currLikelihood)
		currDirection, currLikelihood = checkEast(data, R, P, rStep, pStep, currDirection, currLikelihood)
		currDirection, currLikelihood = checkSouthEast(data, R, P, rStep, pStep, currDirection, currLikelihood)
		currDirection, currLikelihood = checkSouth(data, R, P, rStep, pStep, currDirection, currLikelihood)
		currDirection, currLikelihood = checkSouthWest(data, R, P, rStep, pStep, currDirection, currLikelihood)
	case east:
		currDirection, currLikelihood = checkNorth(data, R, P, rStep, pStep, currDirection, currLikelihood)
		currDirection, currLikelihood = checkNorthEast(data, R, P, rStep, pStep, currDirection, currLikelihood)
		currDirection, currLikelihood = checkEast(data, R, P, rStep, pStep, currDirection, currLikelihood)
		currDirection, currLikelihood = checkSouthEast(data, R, P, rStep, pStep, currDirection, currLikelihood)
		currDirection, currLikelihood = checkSouth(data, R, P, rStep, pStep, currDirection, currLikelihood)
	case northeast:
		currDirection, currLikelihood = checkNorthWest(data, R, P, rStep, pStep, currDirection, currLikelihood)
		currDirection, currLikelihood = checkNorth(data, R, P, rStep, pStep, currDirection, currLikelihood)
		currDirection, currLikelihood = checkNorthEast(data, R, P, rStep, pStep, currDirection, currLikelihood)
		currDirection, currLikelihood = checkEast(data, R, P, rStep, pStep, currDirection, currLikelihood)
		currDirection, currLikelihood = checkSouthEast(data, R, P, rStep, pStep, currDirection, currLikelihood)
	default:
		log.Panic("Error: nextDirection. Unknown direction.\n")
	}

	if currLikelihood == prevLikelihood {
		return neutral, currLikelihood
	}

	return currDirection, currLikelihood
}

// firstDirection determines the first direction for coordinate ascent based on the data and initial parameter set R and P.
// the second return is the likelihood of the data given the initial parameters R and P.
func firstDirection(data []int, R float64, P float64, rStep float64, pStep float64) (direction, float64) {
	var currDirection direction = neutral
	currLikelihood := zeroTruncatedNegativeBinomialLogLikelihood(data, R, P)

	currDirection, currLikelihood = checkNorth(data, R, P, rStep, pStep, currDirection, currLikelihood)
	currDirection, currLikelihood = checkNorthEast(data, R, P, rStep, pStep, currDirection, currLikelihood)
	currDirection, currLikelihood = checkEast(data, R, P, rStep, pStep, currDirection, currLikelihood)
	currDirection, currLikelihood = checkSouthEast(data, R, P, rStep, pStep, currDirection, currLikelihood)
	currDirection, currLikelihood = checkSouth(data, R, P, rStep, pStep, currDirection, currLikelihood)
	currDirection, currLikelihood = checkSouthWest(data, R, P, rStep, pStep, currDirection, currLikelihood)
	currDirection, currLikelihood = checkWest(data, R, P, rStep, pStep, currDirection, currLikelihood)
	currDirection, currLikelihood = checkNorthWest(data, R, P, rStep, pStep, currDirection, currLikelihood)

	return currDirection, currLikelihood
}

// moveInDirection is a helper function of ZeroTruncatedNegativeBinomial. It moves the current parameter set in the
// direction 'currDirection' by one step.
func moveInDirection(R float64, P float64, rStep float64, pStep float64, currDirection direction) (float64, float64) {
	switch currDirection {
	case neutral:
		return R, P
	case north:
		return R + rStep, P
	case northeast:
		if P+pStep < 0.999 {
			return R + rStep, P + pStep
		}
		return R + rStep, P
	case east:
		if P+pStep <= 1 {
			return R, P + pStep
		}
		log.Panic("Error: direct collision with paramter space barrier.\n")
	case southeast:
		if R-rStep < 0.001 && P+pStep > 0.999 {
			log.Panic("Error: direct collision with parameter space barrier.\n")
		}
		if R-rStep < 0.001 {
			return R, P + pStep
		}
		if P+pStep > 0.999 {
			return R - rStep, P
		}
		return R - rStep, P + pStep
	case south:
		if R-rStep < 0.001 {
			log.Panic("Error: direct collision with parameter space barrier.\n")
		}
		return R - rStep, P
	case southwest:
		if R-rStep < 0.001 && P-pStep < 0.001 {
			log.Panic("Error: direct collision with parameter space barrier.\n")
		}
		if R-rStep < 0.001 {
			return R, P - pStep
		}
		return R - rStep, P - pStep
	case west:
		if P-pStep <= 0 {
			log.Panic("Error: direct collision with parameter space barrier.\n")
		}
		return R, P - pStep
	case northwest:
		if P-pStep <= 0 {
			return R + rStep, P
		}
		return R + rStep, P - pStep
	default:
		log.Panic("Error: MoveInDirection: Unknown direction.\n")
		return -1, -1
	}
	return -1, -1
}

// ZeroTruncatedNegativeBinomial fits the maximum likelihood zero-truncated negative binomial distribution to input count data using coordinate ascent.
// Input count data is expected to be in histogram form, where data[i] is equal to the number of observations with a score of i.
// As the fit is to a zero-truncated distribution, any value in data[0] is ignored.
// note that because the loss function is logLikelihood, we want the least negative value, so we are maximizing
// parameter space is defined where R is north/south (south is low, north is high), and P is east/west.
func ZeroTruncatedNegativeBinomial(data []int, R float64, P float64, rStep float64, pStep float64) (float64, float64) {
	if R <= 0 {
		log.Fatalf("Error: initial R value must be greater than 0. Found: %v.\n", R)
	}
	if P <= 0 || P >= 1 {
		log.Fatalf("Error: initial P value must a valid probability. Found: %v.\n", P)
	}

	//fmt.Printf("Getting first direction.\n")
	currDirection, currLikelihood := firstDirection(data, R, P, rStep, pStep)
	//fmt.Printf("Got first direction.\n")

	for currDirection != neutral {
		R, P = moveInDirection(R, P, rStep, pStep, currDirection)
		currDirection, currLikelihood = nextDirection(data, R, P, rStep, pStep, currDirection, currLikelihood)
	}

	return R, P
}

// randNegativeBinomial is a function used only in testing. It generates negative binomial distributed random variates
// from a distribution defined by parameters r and p.
func randNegativeBinomial(r float64, p float64) int {
	s := 0           //successes
	f := 0           //failures
	for s < int(r) { // for fewer than r successes
		if rand.Float64() < p {
			s++
		} else {
			f++
		}
	}
	return f
}
