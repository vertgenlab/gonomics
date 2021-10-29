package phylo

import (
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/numbers"
	"log"
	"math"
)

//AccelDistancesAndWeights represents a set of all observed pairwise distances between the four species used for acceleration analysis.
//While distances are integers, we store them as float64 to avoid casting in other functions as a way to improve readability.
//The six weight fields are the corresponding weights for each distance (w = 1/ (D^2))
type AccelDistancesAndWeights struct {
	DhumChimp   float64
	DhumGor     float64
	DhumOrang   float64
	DchimpGor   float64
	DchimpOrang float64
	DgorOrang   float64
	WhumChimp   float64
	WhumGor     float64
	WhumOrang   float64
	WchimpGor   float64
	WchimpOrang float64
	WgorOrang   float64
}

//AccelBranchLengths represents the set of branch lengths corresponding to a particular distance matrix used for acceleration analysis on a 4-species alingment.
type AccelBranchLengths struct {
	BhumHca   float64
	BchimpHca float64
	BhcaHga   float64
	BhgaGor   float64
	BhgaOrang float64
}

//AccelSubTreeLeft is a tree with three leaves and three branches, joined at the single internal node.
//Dij represent the observed pairwise distances between two species. vi represents the length of the branch between a species i and the internal node at the current stage of optimization.
//Wij are the weight constants corresponding to each distance Dij such that Wij = 1 / (Dij^2)
type AccelSubTreeLeft struct {
	DhumChimp float64
	DhumHga   float64
	DchimpHga float64
	VhumHca   float64
	VchimpHca float64
	VhcaHga   float64
	WhumChimp float64
	WhumHga   float64
	WchimpHga float64
}

//AccelSubTreeRight is a tree with three leaves and three branches, joined at the single internal node.
//Dij represent the observed pairwise distances between two species. vi represents the length of the branch between a species i and the internal node at the current stage of optimization.
//Wij are the weight constants corresponding to each distance Dij such that Wij = 1 / (Dij^2)
type AccelSubTreeRight struct {
	DgorOrang float64
	DhcaGor   float64
	DhcaOrang float64
	VhcaHga   float64
	VhgaGor   float64
	VhgaOrang float64
	WgorOrang float64
	WhcaGor   float64
	WhcaOrang float64
}

//BranchLengthsAlternatingLeastSquares calculates the optimal branch lengths for a given set of distances. See the multiFaAcceleration readme for a detailed theoretical description of this algorithm.
func BranchLengthsAlternatingLeastSquares(d AccelDistancesAndWeights, allowNegative bool, verbose bool, zeroDistanceWeightConstant float64, epsilon float64) AccelBranchLengths {
	var answer = AccelBranchLengths{1, 1, 1, 1, 1}
	var Q float64 = calculateQ(d, answer)
	var nextQ float64
	var currDiff float64 = epsilon + 1 //set currDiff to something larger than epsilon so that we make it into the loop the first time.
	var leftSub AccelSubTreeLeft
	var rightSub AccelSubTreeRight
	var maxIteration, i = 1000, 0
	var oldAnswer = AccelBranchLengths{1, 1, 1, 1, 1}

	for currDiff > epsilon && i < maxIteration {
		oldAnswer = answer
		pruneLeft(d, answer, &leftSub, zeroDistanceWeightConstant)
		optimizeSubtreeLeft(&leftSub, allowNegative, verbose, &answer)
		pruneRight(d, answer, &rightSub, zeroDistanceWeightConstant)
		optimizeSubtreeRight(&rightSub, allowNegative, verbose, &answer)
		nextQ = calculateQ(d, answer)
		//DEBUG: log.Printf("nextQ: %e. currDiff: %e. Here were the branch lengths: %f. %f. %f. %f. %f.", nextQ, currDiff, answer.B1, answer.B2, answer.B3, answer.B4, answer.B5)
		currDiff = math.Abs(Q - nextQ)
		if nextQ > Q { //nextQ is higher than Q, which means we got "worse"
			answer = oldAnswer //we will exit the loop next time, so we want the old answer, which has the lower of the two terminal Q estimates.
			currDiff = 0
		}
		Q = nextQ
		i++
	}
	if i >= maxIteration {
		log.Fatalf("Failed to converge on a tree with these distances. DhumHca: %f, DhumGor: %f, DhumOrang: %f, DchimpGor: %f, DchimpOrang: %f, DgorOrang: %f.", d.DhumChimp, d.DhumGor, d.DhumOrang, d.DchimpGor, d.DchimpOrang, d.DgorOrang)
	}
	return answer
}

//AccelFourWaySnpDistancesAndweights generates an AccelDistancesAndWeights struct from SNP distance, which includes only SNPs.
//Species index order is as follows: 0-Human, 1-Chimp, 2-Gor, 3-Orang
func AccelFourWaySnpDistancesAndWeights(records []fasta.Fasta, alignmentCounter int, windowSize int, d *AccelDistancesAndWeights, zeroDistanceWeightConstant float64) bool {
	//first we clear the values in d.
	d.DhumChimp, d.DhumGor, d.DhumOrang, d.DchimpGor, d.DchimpOrang, d.DgorOrang = 0, 0, 0, 0, 0, 0
	var baseCount, i int = 0, 0
	var reachedEnd bool = false

	if len(records) != 4 {
		log.Fatalf("multiFaAcceleration must take in a four-way multiple alignment.")
	}
	for i = alignmentCounter; baseCount < windowSize && i < len(records[0].Seq); i++ {
		if records[0].Seq[i] != dna.Gap {
			baseCount++
		}
		if IsUngappedColumn(records, i) {
			if records[0].Seq[i] != records[1].Seq[i] {
				d.DhumChimp++
			}
			if records[0].Seq[i] != records[2].Seq[i] {
				d.DhumGor++
			}
			if records[0].Seq[i] != records[3].Seq[i] {
				d.DhumOrang++
			}
			if records[1].Seq[i] != records[2].Seq[i] {
				d.DchimpGor++
			}
			if records[1].Seq[i] != records[3].Seq[i] {
				d.DchimpOrang++
			}
			if records[2].Seq[i] != records[3].Seq[i] {
				d.DgorOrang++
			}
		}
	}
	if baseCount != windowSize {
		reachedEnd = true
	}
	calculateWeights(d, zeroDistanceWeightConstant)
	return reachedEnd
}

//AccelFourWayMutationDistancesAndWeights generates an AccelDistancesAndWeights struct from mutation distances, which includes SNPs and INDELs, where each INDEL counts as one mutation regardless of length.
//Species index order is as follows: 0-Human, 1-Chimp, 2-Gor, 3-Orang
func AccelFourWayMutationDistancesAndWeights(records []fasta.Fasta, alignmentCounter int, windowSize int, d *AccelDistancesAndWeights, zeroDistanceWeightConstant float64) bool {
	//first we clear the values in D.
	d.DhumChimp, d.DhumGor, d.DhumOrang, d.DchimpGor, d.DchimpOrang, d.DgorOrang = 0, 0, 0, 0, 0, 0
	var DhumChimptmp int
	var reachedEnd bool
	var alnEnd int
	DhumChimptmp, reachedEnd, alnEnd = fasta.PairwiseMutationDistanceReferenceWindow(records[0], records[1], alignmentCounter, windowSize)
	d.DhumChimp = float64(DhumChimptmp)
	d.DhumGor = float64(fasta.PairwiseMutationDistanceInRange(records[0], records[2], alignmentCounter, alnEnd))
	d.DhumOrang = float64(fasta.PairwiseMutationDistanceInRange(records[0], records[3], alignmentCounter, alnEnd))
	d.DchimpGor = float64(fasta.PairwiseMutationDistanceInRange(records[1], records[2], alignmentCounter, alnEnd))
	d.DchimpOrang = float64(fasta.PairwiseMutationDistanceInRange(records[1], records[3], alignmentCounter, alnEnd))
	d.DgorOrang = float64(fasta.PairwiseMutationDistanceInRange(records[2], records[3], alignmentCounter, alnEnd))
	calculateWeights(d, zeroDistanceWeightConstant)
	return reachedEnd
}

//calculateWeights is a helper function of the distanceAndWeights functions, produces the weight constant for each distance.
func calculateWeights(d *AccelDistancesAndWeights, zeroDistanceWeightConstant float64) {
	d.WhumChimp = calculateWeight(d.DhumChimp, zeroDistanceWeightConstant)
	d.WhumGor = calculateWeight(d.DhumGor, zeroDistanceWeightConstant)
	d.WhumOrang = calculateWeight(d.DhumOrang, zeroDistanceWeightConstant)
	d.WchimpGor = calculateWeight(d.DchimpGor, zeroDistanceWeightConstant)
	d.WchimpOrang = calculateWeight(d.DchimpOrang, zeroDistanceWeightConstant)
	d.WgorOrang = calculateWeight(d.DgorOrang, zeroDistanceWeightConstant)
}

//calculateWeight is a helper function of calculateWeights, and performs each individual weight calculation.
func calculateWeight(d float64, zeroDistanceWeightConstant float64) float64 {
	if d == 0 {
		return zeroDistanceWeightConstant
	}
	return 1.0 / math.Pow(d, 2)
}

//IsUngappedColumn determines if an alignment column is comprised of bases (not gaps) for each species.
func IsUngappedColumn(records []fasta.Fasta, index int) bool {
	for i := range records {
		if !isUngappedBase(records[i].Seq[index]) {
			return false
		}
	}
	return true
}

//isUngappedBase is a helper function of isUngappedColumn. True if a dna.Base is a base, not an N, gap, or dot.
func isUngappedBase(b dna.Base) bool {
	if b == dna.A || b == dna.T || b == dna.C || b == dna.G {
		return true
	}
	if b == dna.LowerA || b == dna.LowerC || b == dna.LowerG || b == dna.LowerT {
		return true
	}
	return false
}

//A helper function of BranchLengthsAlternatingLeastSquares. Reduce the four species tree to the subtree containing species 0, 1, and the ancestor of 2/3.
func pruneLeft(d AccelDistancesAndWeights, b AccelBranchLengths, sub *AccelSubTreeLeft, ZeroDistanceWeightConstant float64) {
	sub.DhumChimp = d.DhumChimp
	sub.DhumHga = (d.WhumGor*(d.DhumGor-b.BhgaGor) + d.WhumOrang*(d.DhumOrang-b.BhgaOrang)) / (d.WhumGor + d.WhumOrang)
	sub.DchimpHga = (d.WchimpGor*(d.DchimpGor-b.BhgaGor) + d.WchimpOrang*(d.DchimpOrang-b.BhgaOrang)) / (d.WchimpGor + d.WchimpOrang)

	sub.WhumChimp = calculateWeight(sub.DhumChimp, ZeroDistanceWeightConstant)
	sub.WhumHga = calculateWeight(sub.DhumHga, ZeroDistanceWeightConstant)
	sub.WchimpHga = calculateWeight(sub.DchimpHga, ZeroDistanceWeightConstant)
}

//A helper function of BranchLengthsAlternatingLeastSquares. Reduce the four species tree to the subtree containing 2, 3, and the ancestor of 0/1.
func pruneRight(d AccelDistancesAndWeights, b AccelBranchLengths, sub *AccelSubTreeRight, ZeroDistanceWeightConstant float64) {
	sub.DgorOrang = d.DgorOrang
	sub.DhcaGor = (d.WhumGor*(d.DhumGor-b.BhumHca) + d.WchimpGor*(d.DchimpGor-b.BchimpHca)) / (d.WhumGor + d.WchimpGor)
	sub.DhcaOrang = (d.WhumOrang*(d.DhumOrang-b.BhumHca) + d.WchimpOrang*(d.DchimpOrang-b.BchimpHca)) / (d.WhumOrang + d.WchimpOrang)

	sub.WgorOrang = d.DgorOrang
	sub.WhcaGor = sub.DhcaGor
	sub.WhcaOrang = sub.DhcaOrang
}

//optimizeSubtreeLeft is a helper function of BranchLengthsAlternatingLeastSquares and computes the optimal branch lengths for the left subtree at a particular iteration.
func optimizeSubtreeLeft(sub *AccelSubTreeLeft, allowNegative bool, verbose bool, answer *AccelBranchLengths) {
	sub.VhumHca = (sub.DhumChimp + sub.DhumHga - sub.DchimpHga) / 2.0
	sub.VchimpHca = (sub.DhumChimp + sub.DchimpHga - sub.DhumHga) / 2.0
	sub.VhcaHga = (sub.DhumHga + sub.DchimpHga - sub.DhumChimp) / 2.0
	if allowNegative {
		answer.BhumHca = sub.VhumHca
		answer.BchimpHca = sub.VchimpHca
		answer.BhcaHga = sub.VhcaHga
		return
	}
	if sub.VhumHca < 0 && sub.VchimpHca < 0 && sub.VhcaHga < 0 {
		if verbose {
			log.Printf("WARNING: All branches are negative.") //TODO: Should this error out?
		}
		sub.VhumHca = 0
		sub.VchimpHca = 0
		sub.VhcaHga = 0
	} else if sub.VhumHca < 0 {
		sub.VhumHca = 0
		if sub.VchimpHca < 0 {
			sub.VchimpHca = 0
			sub.VhcaHga = nonNegativeApproximation(sub.DhumHga, sub.DchimpHga, sub.VhumHca, sub.VchimpHca, sub.WhumHga, sub.WchimpHga)
		} else if sub.VhcaHga < 0 {
			sub.VhcaHga = 0
			sub.VchimpHca = nonNegativeApproximation(sub.DhumChimp, sub.DchimpHga, sub.VhumHca, sub.VhcaHga, sub.WhumChimp, sub.WchimpHga)
		} else {
			sub.VhcaHga = nonNegativeApproximation(sub.DhumHga, sub.DchimpHga, sub.VhumHca, sub.VchimpHca, sub.WhumHga, sub.WchimpHga)
			sub.VchimpHca = nonNegativeApproximation(sub.DhumChimp, sub.DchimpHga, sub.VhumHca, sub.VhcaHga, sub.WhumChimp, sub.WchimpHga)
		}
	} else if sub.VchimpHca < 0 {
		sub.VchimpHca = 0
		if sub.VhcaHga < 0 {
			sub.VhcaHga = 0
			sub.VhumHca = nonNegativeApproximation(sub.DhumHga, sub.DhumChimp, sub.VhcaHga, sub.VchimpHca, sub.WhumHga, sub.WhumChimp)
		} else {
			sub.VhumHca = nonNegativeApproximation(sub.DhumHga, sub.DhumChimp, sub.VhcaHga, sub.VchimpHca, sub.WhumHga, sub.WhumChimp)
			sub.VhcaHga = nonNegativeApproximation(sub.DhumHga, sub.DchimpHga, sub.VhumHca, sub.VchimpHca, sub.WhumHga, sub.WchimpHga)
		}
	} else if sub.VhcaHga < 0 {
		sub.VhcaHga = 0
		sub.VhumHca = nonNegativeApproximation(sub.DhumHga, sub.DhumChimp, sub.VhcaHga, sub.VchimpHca, sub.WhumHga, sub.WhumChimp)
		sub.VchimpHca = nonNegativeApproximation(sub.DhumChimp, sub.DchimpHga, sub.VhumHca, sub.VhcaHga, sub.WhumChimp, sub.WchimpHga)
	}
	answer.BhumHca = sub.VhumHca
	answer.BchimpHca = sub.VchimpHca
	answer.BhcaHga = sub.VhcaHga
}

//optimizeSubtreeRight is a helper function of BranchLengthsAlternatingLeastSquares and computes the optimal branch lengths for the right subtree at a particular iteration.
func optimizeSubtreeRight(sub *AccelSubTreeRight, allowNegative bool, verbose bool, answer *AccelBranchLengths) {
	sub.VhcaHga = (sub.DhcaGor + sub.DhcaOrang - sub.DgorOrang) / 2.0
	sub.VhgaGor = (sub.DhcaGor + sub.DgorOrang - sub.DhcaOrang) / 2.0
	sub.VhgaOrang = (sub.DhcaOrang + sub.DgorOrang - sub.DhcaGor) / 2.0
	if allowNegative {
		answer.BhcaHga = sub.VhcaHga
		answer.BhgaGor = sub.VhgaGor
		answer.BhgaOrang = sub.VhgaOrang
		return
	}

	if sub.VhcaHga < 0 && sub.VhgaGor < 0 && sub.VhgaOrang < 0 {
		if verbose {
			log.Printf("WARNING: All branches are negative.") //TODO: Should this error out?
		}
		sub.VhcaHga = 0
		sub.VhgaGor = 0
		sub.VhgaOrang = 0
	} else if sub.VhcaHga < 0 {
		sub.VhcaHga = 0
		if sub.VhgaGor < 0 {
			sub.VhgaGor = 0
			sub.VhgaOrang = nonNegativeApproximation(sub.DhcaOrang, sub.DgorOrang, sub.VhcaHga, sub.VhgaGor, sub.WhcaOrang, sub.WgorOrang)
		} else if sub.VhgaOrang < 0 {
			sub.VhgaOrang = 0
			sub.VhgaGor = nonNegativeApproximation(sub.DhcaGor, sub.DgorOrang, sub.VhcaHga, sub.VhgaOrang, sub.WhcaGor, sub.WgorOrang)
		} else {
			sub.VhgaOrang = nonNegativeApproximation(sub.DhcaOrang, sub.DgorOrang, sub.VhcaHga, sub.VhgaGor, sub.WhcaOrang, sub.WgorOrang)
			sub.VhgaGor = nonNegativeApproximation(sub.DhcaGor, sub.DgorOrang, sub.VhcaHga, sub.VhgaOrang, sub.WhcaGor, sub.WgorOrang)
		}
	} else if sub.VhgaGor < 0 {
		sub.VhgaGor = 0
		if sub.VhgaOrang < 0 {
			sub.VhgaOrang = 0
			sub.VhcaHga = nonNegativeApproximation(sub.DhcaGor, sub.DhcaOrang, sub.VhgaGor, sub.VhgaOrang, sub.WhcaGor, sub.WhcaOrang)
		} else {
			sub.VhgaOrang = nonNegativeApproximation(sub.DhcaOrang, sub.DgorOrang, sub.VhcaHga, sub.VhgaGor, sub.WhcaOrang, sub.WgorOrang)
			sub.VhcaHga = nonNegativeApproximation(sub.DhcaGor, sub.DhcaOrang, sub.VhgaGor, sub.VhgaOrang, sub.WhcaGor, sub.WhcaOrang)
		}
	} else if sub.VhgaOrang < 0 {
		sub.VhgaOrang = 0
		sub.VhgaGor = nonNegativeApproximation(sub.DhcaGor, sub.DgorOrang, sub.VhcaHga, sub.VhgaOrang, sub.WhcaGor, sub.WgorOrang)
		sub.VhcaHga = nonNegativeApproximation(sub.DhcaGor, sub.DhcaOrang, sub.VhgaGor, sub.VhgaOrang, sub.WhcaGor, sub.WhcaOrang)
	}
	answer.BhcaHga = sub.VhcaHga
	answer.BhgaGor = sub.VhgaGor
	answer.BhgaOrang = sub.VhgaOrang
}

//a helper function of optimizeSubtree.
//If we constrain branch lengths to be nonNegative, we apply this correction when the minimum Q is achieved at negative branch lengths for a subtree.
func nonNegativeApproximation(d1 float64, d2 float64, v1 float64, v2 float64, w1 float64, w2 float64) float64 {
	return numbers.MaxFloat64((w1*(d1-v1) + w2*(d2-v2)) / (w1 + w2), 0)//ensures the estimate is non-negative
}

//For a set of distances and corresponding branch lengths, determine the value of Q, the Fitch-Margoliash least squares error term.
func calculateQ(d AccelDistancesAndWeights, b AccelBranchLengths) float64 {
	var sum float64 = 0
	sum += d.WhumChimp * math.Pow(d.DhumChimp-(b.BchimpHca+b.BchimpHca), 2)
	sum += d.WhumGor * math.Pow(d.DhumGor-(b.BhumHca+b.BhcaHga+b.BhgaGor), 2)
	sum += d.WhumOrang * math.Pow(d.DhumOrang-(b.BhumHca+b.BhcaHga+b.BhgaOrang), 2)
	sum += d.WchimpGor * math.Pow(d.DchimpGor-(b.BchimpHca+b.BhcaHga+b.BhgaGor), 2)
	sum += d.WchimpOrang * math.Pow(d.DchimpOrang-(b.BchimpHca+b.BhcaHga+b.BhgaOrang), 2)
	sum += d.WgorOrang * math.Pow(d.DgorOrang-(b.BhgaGor+b.BhgaOrang), 2)
	return sum
}
