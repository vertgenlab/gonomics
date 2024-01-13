package sam

import (
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/numbers"
	"github.com/vertgenlab/gonomics/numbers/logspace"
	"github.com/vertgenlab/gonomics/numbers/parse"
	"log"
	"math"
	"strings"
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

// DiploidBaseToBases converts a DiploidBase type into a slice of two dna.Bases.
func DiploidBaseToBases(base DiploidBase) []dna.Base {
	switch base {
	case AA:
		return []dna.Base{dna.A, dna.A}
	case AC:
		return []dna.Base{dna.A, dna.C}
	case AG:
		return []dna.Base{dna.A, dna.G}
	case AT:
		return []dna.Base{dna.A, dna.T}
	case CC:
		return []dna.Base{dna.C, dna.C}
	case CG:
		return []dna.Base{dna.C, dna.G}
	case CT:
		return []dna.Base{dna.C, dna.T}
	case GG:
		return []dna.Base{dna.G, dna.G}
	case GT:
		return []dna.Base{dna.G, dna.T}
	case TT:
		return []dna.Base{dna.T, dna.T}
	case NN:
		return []dna.Base{dna.N, dna.N}
	}
	log.Fatalf("Error: Unrecognized DiploidBase: %v.\n", base)
	return []dna.Base{}
}

// DiploidBaseString formats a DiploidBase type as a string.
func DiploidBaseString(base DiploidBase) string {
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
	log.Fatalf("Error: Unrecognized DiploidBase: %v.", base)
	return ""
}

// DiploidBaseCallFromPile determines which of 10 DiploidBase genotypes
// are present at the position of an input Pile, preconditioned on the assumption of
// a diploid base at that position (in other words, not considering INDELs).
// epsilon represents the error rate, a parameter for the likelihood function.
// lambda represents the rate of cytosine deamination from postmortem hydrolytic degradation.
func DiploidBaseCallFromPile(p Pile, refBase dna.Base, priorCache [][]float64, homozygousCache [][]float64, heterozygousCache [][]float64, ancientCache AncientLikelihoodCache, epsilon float64, lambda float64) DiploidBase {
	if refBase == dna.N { // we will force a return of diploid NN whenever an N is found in the reference.
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
		case dna.N: // should be redundant with catch above
			return NN
		default:
			log.Fatalf("Error: Reference base was not N, A, C, G, or T. Found: %s.\n", dna.BaseToString(refBase))
		}
	}

	var maxDiploid []DiploidBase
	var maxPosterior float64

	switch refBase {
	case dna.A:
		maxDiploid = []DiploidBase{AA}
		if lambda > 0 {
			maxPosterior = logspace.Multiply(ancientBaseLikelihood(aCount, cCount, gCount, tCount, AA, epsilon, lambda, ancientCache), priorCache[dna.A][AA])
		} else {
			maxPosterior = logspace.Multiply(baseLikelihood(aCount, cCount, gCount, tCount, AA, epsilon, homozygousCache, heterozygousCache), priorCache[dna.A][AA])
		}
	case dna.C:
		maxDiploid = []DiploidBase{CC}
		if lambda > 0 {
			maxPosterior = logspace.Multiply(ancientBaseLikelihood(aCount, cCount, gCount, tCount, CC, epsilon, lambda, ancientCache), priorCache[dna.C][CC])
		} else {
			maxPosterior = logspace.Multiply(baseLikelihood(aCount, cCount, gCount, tCount, CC, epsilon, homozygousCache, heterozygousCache), priorCache[dna.C][CC])
		}
	case dna.G:
		maxDiploid = []DiploidBase{GG}
		if lambda > 0 {
			maxPosterior = logspace.Multiply(ancientBaseLikelihood(aCount, cCount, gCount, tCount, GG, epsilon, lambda, ancientCache), priorCache[dna.G][GG])
		} else {
			maxPosterior = logspace.Multiply(baseLikelihood(aCount, cCount, gCount, tCount, GG, epsilon, homozygousCache, heterozygousCache), priorCache[dna.G][GG])
		}
	case dna.T:
		maxDiploid = []DiploidBase{TT}
		if lambda > 0 {
			maxPosterior = logspace.Multiply(ancientBaseLikelihood(aCount, cCount, gCount, tCount, TT, epsilon, lambda, ancientCache), priorCache[dna.T][TT])
		} else {
			maxPosterior = logspace.Multiply(baseLikelihood(aCount, cCount, gCount, tCount, TT, epsilon, homozygousCache, heterozygousCache), priorCache[dna.T][TT])
		}
	default:
		log.Fatalf("Error: Reference base was not N, A, C, G, or T. Found: %s.\n", dna.BaseToString(refBase))
	}
	var geno DiploidBase
	var currPosterior float64
	for geno = 0; geno < 10; geno++ { // for each genotype, encoded as a DiploidBase byte ranging from 0 to 9, inclusive
		if lambda > 0 {
			currPosterior = logspace.Multiply(ancientBaseLikelihood(aCount, cCount, gCount, tCount, geno, epsilon, lambda, ancientCache), priorCache[refBase][geno])
		} else {
			currPosterior = logspace.Multiply(baseLikelihood(aCount, cCount, gCount, tCount, geno, epsilon, homozygousCache, heterozygousCache), priorCache[refBase][geno])
		}
		if currPosterior > maxPosterior {
			maxPosterior = currPosterior
			maxDiploid = maxDiploid[:1] // clear ties
			maxDiploid[0] = geno
		} else if currPosterior == maxPosterior {
			maxDiploid = append(maxDiploid, geno)
		}
	}

	return maxDiploid[numbers.RandIntInRange(0, len(maxDiploid))] // if two genotypes have the same posterior density, pick one at random.
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
	log.Fatalf("Error: Genotype unknown. Found: %v.\n", geno)
	return 0
}

// homozygousLikelihoodExpression is a helper function o baseLikelihood that evaluates the likelihood function for a
// homozygous genotype. CorrectCount is the number of bases observed from the correct base. IncorrectCount is the sum of
// counts from the other three bases.  Epsilon represents the misclassification rate.
// homozygousCache is a [][]float64 that stores the output values for correctCount/incorrectCount pairs already observed.
func homozygousLikelihoodExpression(correctCount int, incorrectCount int, epsilon float64, homozygousCache [][]float64) float64 {
	if correctCount < len(homozygousCache) && incorrectCount < len(homozygousCache[correctCount]) { // if the coverage is within the cache bounds
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
	if correctCount < len(heterozygousCache) && incorrectCount < len(heterozygousCache[correctCount]) { // if the coverage is within the cache bounds
		if heterozygousCache[correctCount][incorrectCount] != 0 { //if this pile has been seen in the cache already
			return heterozygousCache[correctCount][incorrectCount]
		} else {
			s := logspace.Pow(math.Log(0.5-epsilon), float64(correctCount))
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

// MakeDiploidBasePriorCache is a helper function used in samAssembler before running DiploidBaseCallFromPile.
// Constructs a matrix ([][]float64) of form answer[refBase][outGenotype], where index of refBase: 0=a, 1=c, 2=g, 3=t.
// Index of outGenotype follows the indices of the DiploidBase type.
// Answers are log transformed.
func MakeDiploidBasePriorCache(delta float64, gamma float64) [][]float64 {
	Tv := delta / (2.0 + gamma) // the probability of transversions
	Tr := gamma * Tv            // the probability of transitions
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

func MakeDiploidBaseEmpiricalPriorCache(inFile string) [][]float64 {
	lines := fileio.Read(inFile)
	if len(lines) != 5 {
		log.Fatalf("Error: expected five lines in empirical prior file. Found: %v.\n", len(lines))
	}
	wordsA := strings.Split(lines[1], "\t")
	if len(wordsA) != 11 {
		log.Fatalf("Error: expected 11 fields in empirical prior file, row 2. Found: %v.\n", len(wordsA))
	}
	wordsC := strings.Split(lines[2], "\t")
	if len(wordsC) != 11 {
		log.Fatalf("Error: expected 11 fields in empirical prior file, row 3. Found: %v.\n", len(wordsC))
	}
	wordsG := strings.Split(lines[3], "\t")
	if len(wordsG) != 11 {
		log.Fatalf("Error: expected 11 fields in empirical prior file, row 4. Found: %v.\n", len(wordsG))
	}
	wordsT := strings.Split(lines[4], "\t")
	if len(wordsT) != 11 {
		log.Fatalf("Error: expected 11 fields in empirical prior file, row 5. Found: %v.\n", len(wordsT))
	}
	return [][]float64{
		{parse.StringToFloat64(wordsA[1]),
			parse.StringToFloat64(wordsA[2]),
			parse.StringToFloat64(wordsA[3]),
			parse.StringToFloat64(wordsA[4]),
			parse.StringToFloat64(wordsA[5]),
			parse.StringToFloat64(wordsA[6]),
			parse.StringToFloat64(wordsA[7]),
			parse.StringToFloat64(wordsA[8]),
			parse.StringToFloat64(wordsA[9]),
			parse.StringToFloat64(wordsA[10]),
		},
		{parse.StringToFloat64(wordsC[1]),
			parse.StringToFloat64(wordsC[2]),
			parse.StringToFloat64(wordsC[3]),
			parse.StringToFloat64(wordsC[4]),
			parse.StringToFloat64(wordsC[5]),
			parse.StringToFloat64(wordsC[6]),
			parse.StringToFloat64(wordsC[7]),
			parse.StringToFloat64(wordsC[8]),
			parse.StringToFloat64(wordsC[9]),
			parse.StringToFloat64(wordsC[10]),
		},
		{parse.StringToFloat64(wordsG[1]),
			parse.StringToFloat64(wordsG[2]),
			parse.StringToFloat64(wordsG[3]),
			parse.StringToFloat64(wordsG[4]),
			parse.StringToFloat64(wordsG[5]),
			parse.StringToFloat64(wordsG[6]),
			parse.StringToFloat64(wordsG[7]),
			parse.StringToFloat64(wordsG[8]),
			parse.StringToFloat64(wordsG[9]),
			parse.StringToFloat64(wordsG[10]),
		},
		{parse.StringToFloat64(wordsT[1]),
			parse.StringToFloat64(wordsT[2]),
			parse.StringToFloat64(wordsT[3]),
			parse.StringToFloat64(wordsT[4]),
			parse.StringToFloat64(wordsT[5]),
			parse.StringToFloat64(wordsT[6]),
			parse.StringToFloat64(wordsT[7]),
			parse.StringToFloat64(wordsT[8]),
			parse.StringToFloat64(wordsT[9]),
			parse.StringToFloat64(wordsT[10]),
		},
	}
}

// MakeDiploidBaseFlatPriorCache returns an uninformative prior distribution for genotypes.
// With 10 genotypes, each genotype has an equal prior probability of 0.1.
// Answers are log transformed.
func MakeDiploidBaseFlatPriorCache() [][]float64 {
	value := math.Log(0.1)

	return [][]float64{{value, value, value, value, value, value, value, value, value, value},
		{value, value, value, value, value, value, value, value, value, value},
		{value, value, value, value, value, value, value, value, value, value},
		{value, value, value, value, value, value, value, value, value, value}}
}
