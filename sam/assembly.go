package sam

import (
	"fmt"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/numbers"
	"github.com/vertgenlab/gonomics/numbers/logspace"
	"log"
	"math"
)

type DiploidBase byte

const (
	AA DiploidBase = 0
	AC DiploidBase = 1
	AG DiploidBase = 2
	AT DiploidBase = 3
	CC DiploidBase = 4
	CG DiploidBase = 5
	CT DiploidBase = 6
	GG DiploidBase = 7
	GT DiploidBase = 8
	TT DiploidBase = 9
	NN DiploidBase = 10
)

// diploidBaseString formats a DiploidBase struct as a string for debugging.
func diploidBaseString(base DiploidBase) string {
	switch base {
	case AA:
		return "AA"
	case AC:
		return "AC"
	case AG:
		return "AG"
	case AT:
		return "AT"
	case CC:
		return "CC"
	case CG:
		return "CG"
	case CT:
		return "CT"
	case GG:
		return "GG"
	case GT:
		return "GT"
	case TT:
		return "TT"
	case NN:
		return "NN"
	}
	log.Fatalf("Unrecognized DiploidBase: %v.", base)
	return ""
}

//DiploidBaseCallFromPile determines which of 10 DiploidBase genotypes
// are present at the position of an input Pile, preconditioned on the assumption of
// a diploid base at that position (in other words, not considering INDELs).
// epsilon represents the error rate, a parameter for the likelihood function.
func DiploidBaseCallFromPile(p Pile, refBase dna.Base, priorCache [][]float64, homozygousCache [][]float64, heterozygousCache [][]float64, epsilon float64) DiploidBase {
	if refBase == dna.N { //we will force a return of diploid NN whenever an N is found in the reference.
		return NN
	}
	var aCount = p.CountF[dna.A] + p.CountR[dna.A]
	var cCount = p.CountF[dna.C] + p.CountR[dna.C]
	var gCount = p.CountF[dna.G] + p.CountR[dna.G]
	var tCount = p.CountF[dna.T] + p.CountR[dna.T]
	var baseCount = aCount + cCount + gCount + tCount

	if baseCount < 1 {
		switch refBase {
		case dna.A:
			return AA
		case dna.C:
			return CC
		case dna.G:
			return GG
		case dna.T:
			return TT
		case dna.N: //should be redundant with catch above
			return NN
		default:
			log.Fatalf("Reference base was not N, A, C, G, or T. Found: %s.\n", dna.BaseToString(refBase))
		}
	}

	var maxDiploid []DiploidBase
	var maxPosterior float64

	switch refBase {
	case dna.A:
		maxDiploid = []DiploidBase{AA}
		maxPosterior = logspace.Multiply(baseLikelihood(aCount, cCount, gCount, tCount, AA, epsilon, homozygousCache, heterozygousCache), priorCache[dna.A][AA])
	case dna.C:
		maxDiploid = []DiploidBase{CC}
		maxPosterior = logspace.Multiply(baseLikelihood(aCount, cCount, gCount, tCount, CC, epsilon, homozygousCache, heterozygousCache), priorCache[dna.C][CC])
	case dna.G:
		maxDiploid = []DiploidBase{GG}
		maxPosterior = logspace.Multiply(baseLikelihood(aCount, cCount, gCount, tCount, GG, epsilon, homozygousCache, heterozygousCache), priorCache[dna.G][GG])
	case dna.T:
		maxDiploid = []DiploidBase{TT}
		maxPosterior = logspace.Multiply(baseLikelihood(aCount, cCount, gCount, tCount, TT, epsilon, homozygousCache, heterozygousCache), priorCache[dna.T][TT])
	default:
		log.Fatalf("Reference base was not N, A, C, G, or T. Found: %s.\n", dna.BaseToString(refBase))
	}

	//DEBUG: fmt.Printf("HomoRefBase: %s. Posterior: %v.\n", diploidBaseString(maxDiploid[0]), maxPosterior)

	var geno DiploidBase
	var currPosterior float64
	for geno = 0; geno < 10; geno++ { //for each genotype, encoded as a DiploidBase byte ranging from 0 to 9, inclusive
		currPosterior = logspace.Multiply(baseLikelihood(aCount, cCount, gCount, tCount, geno, epsilon, homozygousCache, heterozygousCache), priorCache[refBase][geno])
		fmt.Printf("CurrGeno: %s. CurrPosterior: %v.\n", diploidBaseString(geno), currPosterior)
		if currPosterior > maxPosterior {
			maxPosterior = currPosterior
			maxDiploid = maxDiploid[:1] //clear ties
			maxDiploid[0] = geno
		} else if currPosterior == maxPosterior {
			maxDiploid = append(maxDiploid, geno)
		}
	}

	return maxDiploid[numbers.RandIntInRange(0, len(maxDiploid))] //if two genotypes have the same posterior density, pick one at random.
}

// baseLikelihood is a helper function of DiploidBaseCallFromPile. For a given genotype (geno), and given counts for the
// four bases, and a given epsilon (misclassification rate), it calculates the genotype likelihood.
// The two caches store values for homozygous and heterozygous genotype pile counts already observed to savce on computation.
func baseLikelihood(aCount int, cCount int, gCount int, tCount int, geno DiploidBase, epsilon float64, homozygousCache [][]float64, heterozygousCache [][]float64) float64 {
	switch geno {
	case AA:
		return homozygousLikelihoodExpression(aCount, cCount+gCount+tCount, epsilon, homozygousCache)
	case AC:
		return heterozygousLikelihoodExpression(aCount+cCount, gCount+tCount, epsilon, heterozygousCache)
	case AG:
		return heterozygousLikelihoodExpression(aCount+gCount, cCount+tCount, epsilon, heterozygousCache)
	case AT:
		return heterozygousLikelihoodExpression(aCount+tCount, cCount+gCount, epsilon, heterozygousCache)
	case CC:
		return homozygousLikelihoodExpression(cCount, aCount+gCount+tCount, epsilon, homozygousCache)
	case CG:
		return heterozygousLikelihoodExpression(cCount+gCount, aCount+tCount, epsilon, heterozygousCache)
	case CT:
		return heterozygousLikelihoodExpression(cCount+tCount, aCount+gCount, epsilon, heterozygousCache)
	case GG:
		return homozygousLikelihoodExpression(gCount, aCount+cCount+tCount, epsilon, homozygousCache)
	case GT:
		return heterozygousLikelihoodExpression(gCount+tCount, aCount+cCount, epsilon, heterozygousCache)
	case TT:
		return homozygousLikelihoodExpression(tCount, aCount+cCount+gCount, epsilon, homozygousCache)
	}
	log.Fatalf("Genotype unknown. Found: %v.\n", geno)
	return 0
}

// homozygousLikelihoodExpression is a helper function o baseLikelihood that evaluates the likelihood function for a
// homozygous genotype. CorrectCount is the number of bases observed from the correct base. IncorrectCount is the sum of
// counts from the other three bases.  Epsilon represents the misclassification rate.
// homozygousCache is a [][]float64 that stores the output values for correctCount/incorrectCount pairs already observed.
func homozygousLikelihoodExpression(correctCount int, incorrectCount int, epsilon float64, homozygousCache [][]float64) float64 {
	if correctCount < len(homozygousCache) && incorrectCount < len(homozygousCache[correctCount]) { //if the coverage is within the cache bounds
		if homozygousCache[correctCount][incorrectCount] != 0 { //if this pile has been seen in the cache already
			return homozygousCache[correctCount][incorrectCount]
		} else {
			s := logspace.Pow(math.Log(1.0-epsilon), float64(correctCount))
			f := logspace.Pow(math.Log(epsilon/3.0), float64(incorrectCount))
			homozygousCache[correctCount][incorrectCount] = logspace.Multiply(s, f)
			return homozygousCache[correctCount][incorrectCount]
		}
	} else { //if the coverage his higher than the cache bounds, calculate by hand
		s := logspace.Pow(math.Log(1.0-epsilon), float64(correctCount))
		f := logspace.Pow(math.Log(epsilon/3.0), float64(incorrectCount))
		return logspace.Multiply(s, f)
	}
}

// heterozygousLikelihoodExpression is a helper function of baseLikelihood that evaluates the likelihood function for a
// heterozygous genotype. CorrectCount is the number of bases observed from the correct two bases.
// IncorrectCount is the count for the other two bases. Epsilon represents the misclassification rate.
// heterozygousCache is a [][]float64 that stores the output values for correctCount/incorrectCount pairs already observed.
func heterozygousLikelihoodExpression(correctCount int, incorrectCount int, epsilon float64, heterozygousCache [][]float64) float64 {
	if correctCount < len(heterozygousCache) && incorrectCount < len(heterozygousCache[correctCount]) { //if the coverage is within the cache bounds
		if heterozygousCache[correctCount][incorrectCount] != 0 { //if this pile has been seen in the cache already
			return heterozygousCache[correctCount][incorrectCount]
		} else {
			s := logspace.Pow(math.Log(1.0-epsilon), float64(correctCount))
			f := logspace.Pow(math.Log(epsilon/3.0), float64(incorrectCount))
			heterozygousCache[correctCount][incorrectCount] = logspace.Multiply(s, f)
			return heterozygousCache[correctCount][incorrectCount]
		}
	} else {
		s := logspace.Pow(math.Log(0.5-(epsilon/3.0)), float64(correctCount))
		f := logspace.Pow(math.Log(epsilon/3.0), float64(incorrectCount))
		return logspace.Multiply(s, f)
	}
}

// makePriorCache is a helper function used in samAssembler before running DiploidBaseCallFromPile.
// Constructs a matrix ([][]float64) of form answer[refBase][outGenotype], where index of refBase: 0=a, 1=c, 2=g, 3=t.
// Index of outGenotype follows the indices of the DiploidBase type.
// Answers are log transformed.
func makePriorCache(delta float64, gamma float64) [][]float64 {
	Tv := delta / (2.0 + gamma) //the probability of transversions
	Tr := gamma * Tv            //the probability of transitions
	oneMinusDeltaSquared := math.Log(math.Pow(1-delta, 2))
	tvSquared := math.Log(Tv * Tv)
	trSquared := math.Log(Tr * Tr)
	twoTvTr := math.Log(Tv * Tr)
	twoTvSquared := math.Log(2 * Tv * Tv)
	twoOneMinusDeltaTv := math.Log(2 * (1 - delta) * Tv)
	twoOneMinusDeltaTr := math.Log(2 * (1 - delta) * Tr)

	return [][]float64{{oneMinusDeltaSquared, twoOneMinusDeltaTv, twoOneMinusDeltaTr, twoOneMinusDeltaTv, tvSquared, twoTvTr, twoTvSquared, trSquared, twoTvTr, tvSquared},
		{tvSquared, twoOneMinusDeltaTv, twoTvSquared, twoTvTr, oneMinusDeltaSquared, twoOneMinusDeltaTv, twoOneMinusDeltaTr, tvSquared, twoTvTr, trSquared},
		{trSquared, twoTvTr, twoOneMinusDeltaTr, twoTvTr, tvSquared, twoOneMinusDeltaTv, twoTvSquared, oneMinusDeltaSquared, twoOneMinusDeltaTv, tvSquared},
		{tvSquared, twoTvTr, twoTvSquared, twoOneMinusDeltaTv, trSquared, twoTvTr, twoOneMinusDeltaTr, twoTvSquared, twoOneMinusDeltaTv, oneMinusDeltaSquared}}
}
