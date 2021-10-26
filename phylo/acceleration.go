package phylo

import (
	"log"
	"math"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/numbers"
)

//AccelDistances represents a set of all observed pairwise distances between the four species used for acceleration analysis.
//While these numbers must be integers, we store them as float64 to avoid casting in other functions as a way to improve readability.
type AccelDistances struct {
	D01 float64
	D02 float64
	D03 float64
	D12 float64
	D13 float64
	D23 float64
}

//AccelBranchLengths represents the set of branch lengths corresponding to a particular distance matrix used for acceleration analysis on a 4-species alingment.
type AccelBranchLengths struct {
	B1 float64
	B2 float64
	B3 float64
	B4 float64
	B5 float64
}

//SubTree is a tree with three leaves and three branches, joined at the single internal node.
//Dij represent the observed pairwise distances between two species. vi represents the length of the branch between a species i and the internal node at the current stage of optimization.
type AccelSubTree struct {
	Dab float64
	Dac float64
	Dbc float64
	Va  float64
	Vb  float64
	Vc  float64
}

//BranchLengthsAlternatingLeastSquares calculates the optimal branch lengths for a given set of distances. See the multiFaAcceleration readme for a detailed theoretical description of this algorithm.
func BranchLengthsAlternatingLeastSquares(d AccelDistances, allowNegative bool, verbose bool, zeroDistanceWeightConstant float64, epsilon float64) AccelBranchLengths {
	var answer = AccelBranchLengths{1, 1, 1, 1, 1}
	var Q float64 = calculateQ(d, answer, zeroDistanceWeightConstant)
	var nextQ float64
	var currDiff float64 = epsilon + 1 //set currDiff to something larger than epsilon so that we make it into the loop the first time.
	var sub AccelSubTree
	var maxIteration, i = 1000, 0
	var oldAnswer = AccelBranchLengths{1, 1, 1, 1, 1}

	for currDiff > epsilon && i < maxIteration {
		oldAnswer = answer
		pruneLeft(d, answer, &sub, zeroDistanceWeightConstant)
		answer.B1, answer.B2, answer.B3 = optimizeSubtree(&sub, allowNegative, verbose, zeroDistanceWeightConstant)
		pruneRight(d, answer, &sub, zeroDistanceWeightConstant)
		answer.B4, answer.B5, answer.B3 = optimizeSubtree(&sub, allowNegative, verbose, zeroDistanceWeightConstant)
		nextQ = calculateQ(d, answer, zeroDistanceWeightConstant)
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
		log.Fatalf("Failed to converge on a tree with these distances. D01: %f, D02: %f, D03: %f, D12: %f, D13: %f, D23: %f.", d.D01, d.D02, d.D03, d.D12, d.D13, d.D23)
	}
	return answer
}

//AccelFourWaySnpDistances generates an AccelDistances struct from SNP distance, which includes only SNPs.
func AccelFourWaySnpDistances(records []fasta.Fasta, alignmentCounter int, windowSize int, d *AccelDistances) bool {
	//first we clear the values in d.
	d.D01, d.D02, d.D03, d.D12, d.D13, d.D23 = 0, 0, 0, 0, 0, 0
	var baseCount, i int = 0, 0
	var reachedEnd bool = false

	if len(records) != 4 {
		log.Fatalf("multiFaAcceleration must take in a four-way multiple alignment.")
	}
	for i = alignmentCounter; baseCount < windowSize && i < len(records[0].Seq); i++ {
		if records[0].Seq[i] != dna.Gap {
			baseCount++
		}
		if isUngappedColumn(records, i) {
			if records[0].Seq[i] != records[1].Seq[i] {
				d.D01++
			}
			if records[0].Seq[i] != records[2].Seq[i] {
				d.D02++
			}
			if records[0].Seq[i] != records[3].Seq[i] {
				d.D03++
			}
			if records[1].Seq[i] != records[2].Seq[i] {
				d.D12++
			}
			if records[1].Seq[i] != records[3].Seq[i] {
				d.D13++
			}
			if records[2].Seq[i] != records[3].Seq[i] {
				d.D23++
			}
		}
	}
	if baseCount != windowSize {
		reachedEnd = true
	}
	return reachedEnd
}


//AccelFourWayMutationDistances generates an AccelDistances struct from mutation distances, which includes SNPs and INDELs, where each INDEL counts as one mutation regardless of length.
func AccelFourWayMutationDistances(records []fasta.Fasta, alignmentCounter int, windowSize int, D *AccelDistances) bool {
	//first we clear the values in D.
	D.D01, D.D02, D.D03, D.D12, D.D13, D.D23 = 0, 0, 0, 0, 0, 0
	var D01tmp int
	var reachedEnd bool
	var alnEnd int
	D01tmp, reachedEnd, alnEnd = fasta.PairwiseMutationDistanceReferenceWindow(records[0], records[1], alignmentCounter, windowSize)
	D.D01 = float64(D01tmp)
	D.D02 = float64(fasta.PairwiseMutationDistanceInRange(records[0], records[2], alignmentCounter, alnEnd))
	D.D03 = float64(fasta.PairwiseMutationDistanceInRange(records[0], records[3], alignmentCounter, alnEnd))
	D.D12 = float64(fasta.PairwiseMutationDistanceInRange(records[1], records[2], alignmentCounter, alnEnd))
	D.D13 = float64(fasta.PairwiseMutationDistanceInRange(records[1], records[3], alignmentCounter, alnEnd))
	D.D23 = float64(fasta.PairwiseMutationDistanceInRange(records[2], records[3], alignmentCounter, alnEnd))
	return reachedEnd
}

//a helper function of fourWaySnpDistances, determines if an alignment column is comprised of bases (not gaps) for each species.
func isUngappedColumn(records []fasta.Fasta, index int) bool {
	for i := range records {
		if !isUngappedBase(records[i].Seq[index]) {
			return false
		}
	}
	return true
}

//a helper function of isUngappedColumn. True if a dna.Base is a base, not an N, gap, or dot.
func isUngappedBase(b dna.Base) bool {
	if b == dna.A || b == dna.T || b == dna.C || b == dna.G {
		return true
	}
	if b == dna.LowerA || b == dna.LowerC || b == dna.LowerG || b == dna.LowerT {
		return true
	}
	return false
}


//a helper function of BranchLengthsAlternatingLeastSquares. Calculates the optimal branch lengths for the three branches in a subtree.
func optimizeSubtree(sub *AccelSubTree, allowNegative bool, verbose bool, zeroDistanceWeightConstant float64) (float64, float64, float64) {
	sub.Va = (sub.Dab + sub.Dac - sub.Dbc) / 2.0
	sub.Vb = (sub.Dab + sub.Dbc - sub.Dac) / 2.0
	sub.Vc = (sub.Dac + sub.Dbc - sub.Dac) / 2.0

	if allowNegative {
		return sub.Va, sub.Vb, sub.Vc
	}
	if sub.Va < 0 && sub.Vb < 0 && sub.Vc < 0 {
		if verbose {
			log.Printf("WARNING: All branches are negative.") //TODO: Should this error out?
		}
		sub.Va, sub.Vb, sub.Vc = 0, 0, 0
	} else if sub.Va < 0 && sub.Vb < 0 {
		sub.Va = 0
		sub.Vb = 0
		sub.Vc = nonNegativeApproximation(sub.Dac, sub.Dbc, sub.Va, sub.Vb, zeroDistanceWeightConstant)
	} else if sub.Va < 0 && sub.Vc < 0 {
		sub.Va = 0
		sub.Vc = 0
		sub.Vb = nonNegativeApproximation(sub.Dbc, sub.Dab, sub.Vc, sub.Va, zeroDistanceWeightConstant)
	} else if sub.Vb < 0 && sub.Vc < 0 {
		sub.Vb = 0
		sub.Vc = 0
		sub.Va = nonNegativeApproximation(sub.Dab, sub.Dac, sub.Vb, sub.Vc, zeroDistanceWeightConstant)
	} else if sub.Va < 0 {
		sub.Va = 0
		sub.Vb = nonNegativeApproximation(sub.Dab, sub.Dbc, sub.Va, sub.Vc, zeroDistanceWeightConstant)
		sub.Vc = nonNegativeApproximation(sub.Dac, sub.Dbc, sub.Va, sub.Vb, zeroDistanceWeightConstant)
	} else if sub.Vb < 0 {
		sub.Vb = 0
		sub.Va = nonNegativeApproximation(sub.Dab, sub.Dac, sub.Vb, sub.Vc, zeroDistanceWeightConstant)
		sub.Vc = nonNegativeApproximation(sub.Dab, sub.Dbc, sub.Va, sub.Vc, zeroDistanceWeightConstant)
	} else if sub.Vc < 0 {
		sub.Vc = 0
		sub.Va = nonNegativeApproximation(sub.Dab, sub.Dac, sub.Vb, sub.Vc, zeroDistanceWeightConstant)
		sub.Vb = nonNegativeApproximation(sub.Dab, sub.Dbc, sub.Va, sub.Vc, zeroDistanceWeightConstant)
	}
	return sub.Va, sub.Vb, sub.Vc
}

//a helper function of optimizeSubtree.
//If we constrain branch lengths to be nonNegative, we apply this correction when the minimum Q is achieved at negative branch lengths for a subtree.
func nonNegativeApproximation(d1 float64, d2 float64, v1 float64, v2 float64, ZeroDistanceWeightConstant float64) float64 {
	if d1 == 0 {
		if d2 == 0 {
			return numbers.MaxFloat64(0, (ZeroDistanceWeightConstant*(d1-v1)+ZeroDistanceWeightConstant*(d2-v2))/(2*ZeroDistanceWeightConstant))
		} else {
			return numbers.MaxFloat64(0, ZeroDistanceWeightConstant*(d1-v1)+(1.0/math.Pow(d2, 2)*(d2-v2))) / (ZeroDistanceWeightConstant + (1.0 / math.Pow(d2, 2)))
		}
	} else if d2 == 0 {
		return numbers.MaxFloat64(0, (1.0/(math.Pow(d1, 2))*(d1-v1)+ZeroDistanceWeightConstant*(d2-v2))/((1.0/math.Pow(d1, 2))+ZeroDistanceWeightConstant))
	}
	return numbers.MaxFloat64(0, (1.0/(math.Pow(d1, 2))*(d1-v1)+(1.0/math.Pow(d2, 2)*(d2-v2)))/((1.0/math.Pow(d1, 2))+(1.0/math.Pow(d2, 2))))
}

//A helper function of BranchLengthsAlternatingLeastSquares. Reduce the four species tree to the subtree containing species 0, 1, and the ancestor of 2/3.
func pruneLeft(d AccelDistances, b AccelBranchLengths, sub *AccelSubTree, ZeroDistanceWeightConstant float64) {
	sub.Dab = d.D01
	if d.D03 == 0 {
		if d.D02 == 0 {
			sub.Dac = (ZeroDistanceWeightConstant*(d.D02-b.B4) + ZeroDistanceWeightConstant*(d.D03-b.B5)) / (2 * ZeroDistanceWeightConstant)
		} else {
			sub.Dac = ((1.0/math.Pow(d.D02, 2))*(d.D02-b.B4) + (ZeroDistanceWeightConstant * (d.D03 - b.B5))) / (ZeroDistanceWeightConstant + (1.0 / math.Pow(d.D02, 2)))
		}
	} else if d.D02 == 0 {
		sub.Dac = (ZeroDistanceWeightConstant*(d.D02-b.B4) + ((1.0 / math.Pow(d.D03, 2)) * (d.D03 - b.B5))) / ((1.0 / math.Pow(d.D03, 2)) + ZeroDistanceWeightConstant)
	} else {
		sub.Dac = ((1.0/math.Pow(d.D02, 2))*(d.D02-b.B4) + (1.0/math.Pow(d.D03, 2))*(d.D03-b.B5)) / ((1.0 / math.Pow(d.D03, 2)) + (1.0 / math.Pow(d.D02, 2)))
	}

	if d.D13 == 0 {
		if d.D12 == 0 {
			sub.Dbc = (ZeroDistanceWeightConstant*(d.D12-b.B4) + ZeroDistanceWeightConstant*(d.D13-b.B5)) / (2 * ZeroDistanceWeightConstant)
		} else {
			sub.Dbc = (ZeroDistanceWeightConstant*(d.D13-b.B5) + (1.0/math.Pow(d.D12, 2))*(d.D12-b.B4)) / ((1.0 / math.Pow(d.D12, 2)) + ZeroDistanceWeightConstant)
		}
	} else if d.D12 == 0 {
		sub.Dbc = (ZeroDistanceWeightConstant*(d.D12-b.B4) + (1.0/math.Pow(d.D13, 2))*(d.D13-b.B5)) / ((1.0 / math.Pow(d.D13, 2)) + ZeroDistanceWeightConstant)
	} else {
		sub.Dbc = ((1.0/math.Pow(d.D12, 2))*(d.D12-b.B4) + (1.0/math.Pow(d.D13, 2))*(d.D13-b.B5)) / ((1.0 / math.Pow(d.D13, 2)) + (1.0 / math.Pow(d.D12, 2)))
	}
}

//A helper function of BranchLengthsAlternatingLeastSquares. Reduce the four species tree to the subtree containing 2, 3, and the ancestor of 0/1.
func pruneRight(d AccelDistances, b AccelBranchLengths, sub *AccelSubTree, ZeroDistanceWeightConstant float64) {
	sub.Dac = d.D23

	if d.D02 == 0 {
		if d.D12 == 0 {
			sub.Dac = (ZeroDistanceWeightConstant*(d.D02-b.B1) + ZeroDistanceWeightConstant*(d.D12-b.B2)) / (2.0 * ZeroDistanceWeightConstant)
		} else {
			sub.Dac = ((1.0/math.Pow(d.D12, 2))*(d.D12-b.B2) + ZeroDistanceWeightConstant*(d.D02-b.B1)) / ((1.0 / math.Pow(d.D12, 2)) + ZeroDistanceWeightConstant)
		}
	} else if d.D12 == 0 {
		sub.Dac = ((1.0/math.Pow(d.D02, 2))*(d.D02-b.B1) + ZeroDistanceWeightConstant*(d.D12-b.B2)) / ((1.0 / math.Pow(d.D02, 2)) + ZeroDistanceWeightConstant)
	} else {
		sub.Dab = ((1.0/math.Pow(d.D02, 2))*(d.D02-b.B1) + (1.0/math.Pow(d.D12, 2))*(d.D12-b.B2)) / ((1.0 / math.Pow(d.D02, 2)) + (1.0 / math.Pow(d.D12, 2)))
	}

	if d.D03 == 0 {
		if d.D13 == 0 {
			sub.Dbc = (ZeroDistanceWeightConstant*(d.D03-b.B1) + ZeroDistanceWeightConstant*(d.D13-b.B2)) / (2.0 * ZeroDistanceWeightConstant)
		} else {
			sub.Dbc = ((1.0/math.Pow(d.D13, 2))*(d.D13-b.B2) + ZeroDistanceWeightConstant*(d.D03-b.B1)) / ((1.0 / math.Pow(d.D13, 2)) + ZeroDistanceWeightConstant)
		}
	} else if d.D13 == 0 {
		sub.Dbc = ((1.0/math.Pow(d.D03, 2))*(d.D03-b.B1) + ZeroDistanceWeightConstant*(d.D13-b.B2)) / ((1.0 / math.Pow(d.D03, 2)) + ZeroDistanceWeightConstant)
	} else {
		sub.Dbc = ((1.0/math.Pow(d.D03, 2))*(d.D03-b.B1) + (1.0/math.Pow(d.D13, 2))*(d.D13-b.B2)) / ((1.0 / math.Pow(d.D03, 2)) + (1.0 / math.Pow(d.D13, 2)))
	}
}

//For a set of distances and corresponding branch lengths, determine the value of Q, the Fitch-Margoliash least squares error term.
func calculateQ(d AccelDistances, b AccelBranchLengths, zeroDistanceWeightConstant float64) float64 {
	var sum float64 = 0
	if d.D01 != 0 { //avoid divide by zero error
		sum += math.Pow(d.D01-b.B1-b.B2, 2) / math.Pow(d.D01, 2)
	} else {
		sum += math.Pow(d.D01-b.B1-b.B2, 2) * zeroDistanceWeightConstant
	}
	if d.D02 != 0 {
		sum += math.Pow(d.D02-b.B1-b.B3-b.B4, 2) / math.Pow(d.D02, 2)
	} else {
		sum += math.Pow(d.D02-b.B1-b.B3-b.B4, 2) * zeroDistanceWeightConstant
	}
	if d.D03 != 0 {
		sum += math.Pow(d.D03-b.B1-b.B3-b.B5, 2) / math.Pow(d.D03, 2)
	} else {
		sum += math.Pow(d.D03-b.B1-b.B3-b.B5, 2) * zeroDistanceWeightConstant
	}
	if d.D12 != 0 {
		sum += math.Pow(d.D12-b.B2-b.B3-b.B4, 2) / math.Pow(d.D12, 2)
	} else {
		sum += math.Pow(d.D12-b.B2-b.B3-b.B4, 2) * zeroDistanceWeightConstant
	}
	if d.D13 != 0 {
		sum += math.Pow(d.D13-b.B2-b.B3-b.B5, 2) / math.Pow(d.D13, 2)
	} else {
		sum += math.Pow(d.D13-b.B2-b.B3-b.B5, 2) * zeroDistanceWeightConstant
	}
	if d.D23 != 0 {
		sum += math.Pow(d.D23-b.B4-b.B5, 2) / math.Pow(d.D23, 2)
	} else {
		sum += math.Pow(d.D23-b.B4-b.B5, 2) * zeroDistanceWeightConstant
	}
	return sum
}


