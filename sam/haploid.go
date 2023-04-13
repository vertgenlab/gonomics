package sam

import (
	"math"

	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/numbers"
	"github.com/vertgenlab/gonomics/numbers/logspace"
)

// HaploidCall is a struct to represent the data of a haploid genotype called
// from a Pileup from a haploid genome.
type HaploidCall struct {
	Base      dna.Base
	Insertion string
	Deletion  int
}

// HaploidCallFromPile calls a haploid variant (including Base identity at the Pile position, and INDELs immediately after the position),
// for an input Pile.
func HaploidCallFromPile(p Pile, refBase dna.Base, epsilon float64, haploidBasePriorCache [][]float64, haploidIndelPriorCache []float64, homoBaseCache [][]float64, heteroBaseCache [][]float64, homoIndelCache [][]float64) HaploidCall {
	var answer HaploidCall = HaploidCall{Base: refBase, Insertion: "", Deletion: 0}
	var aCount = p.CountF[dna.A] + p.CountR[dna.A]
	var cCount = p.CountF[dna.C] + p.CountR[dna.C]
	var gCount = p.CountF[dna.G] + p.CountR[dna.G]
	var tCount = p.CountF[dna.T] + p.CountR[dna.T]
	var N = aCount + cCount + gCount + tCount

	if refBase != dna.N { //if we have a real base, we'll calculate posteriors.
		var maxBase = []dna.Base{dna.A}
		var maxPosterior = logspace.Multiply(baseLikelihood(aCount, cCount, gCount, tCount, AA, epsilon, homoBaseCache, heteroBaseCache), haploidBasePriorCache[refBase][dna.A])
		var currPosterior float64

		//CC
		currPosterior = logspace.Multiply(baseLikelihood(aCount, cCount, gCount, tCount, CC, epsilon, homoBaseCache, heteroBaseCache), haploidBasePriorCache[refBase][dna.C])
		if currPosterior > maxPosterior {
			maxBase = maxBase[:1] //clear ties
			maxBase[0] = dna.C
			maxPosterior = currPosterior
		} else if currPosterior == maxPosterior {
			maxBase = append(maxBase, dna.C)
		}
		//GG
		currPosterior = logspace.Multiply(baseLikelihood(aCount, cCount, gCount, tCount, GG, epsilon, homoBaseCache, heteroBaseCache), haploidBasePriorCache[refBase][dna.G])
		if currPosterior > maxPosterior {
			maxBase = maxBase[:1] //clear ties
			maxBase[0] = dna.G
			maxPosterior = currPosterior
		} else if currPosterior == maxPosterior {
			maxBase = append(maxBase, dna.G)
		}
		//TT
		currPosterior = logspace.Multiply(baseLikelihood(aCount, cCount, gCount, tCount, TT, epsilon, homoBaseCache, heteroBaseCache), haploidBasePriorCache[refBase][dna.T])
		if currPosterior > maxPosterior {
			maxBase = maxBase[:1] //clear ties
			maxBase[0] = dna.T
			maxPosterior = currPosterior
		} else if currPosterior == maxPosterior {
			maxBase = append(maxBase, dna.T)
		}
		answer.Base = maxBase[numbers.RandIntInRange(0, len(maxBase))]
	}

	//next we evaluate haploid insertion. We only consider Ia
	var iTot int
	var insertions = make(map[string]int, 0) //this will be the combined insertion depth map, merging the for and rev observations
	for key, value := range p.InsCountF {
		iTot += value
		insertions[key] = value
	}
	for key, value := range p.InsCountR {
		iTot += value
		if _, ok := insertions[key]; ok { // if insertion on reverse strand was seen already on positive strand
			insertions[key] = insertions[key] + value //add reverse read counts to forward read counts
		} else {
			insertions[key] = value
		}
	}
	var IaKey string
	var IaValue int
	for key, value := range insertions {
		if value > IaValue {
			IaKey = key
			IaValue = value
		}
	}
	var B = N - iTot
	var noInsertionPosterior float64
	var insertionPosterior float64
	if IaValue > 0 { //if IaValue == 0, keep initial answer insertion value of ""
		noInsertionPosterior = logspace.Multiply(homozygousIndelLikelihoodExpression(B, IaValue, epsilon, homoIndelCache), haploidIndelPriorCache[0])
		insertionPosterior = logspace.Multiply(homozygousIndelLikelihoodExpression(IaValue, B, epsilon, homoIndelCache), haploidIndelPriorCache[1])
		if insertionPosterior > noInsertionPosterior {
			answer.Insertion = IaKey
		}
	}

	//finally we evaluate haploid deletion. We only consider Da
	var dTot int
	var deletions = make(map[int]int, 0) //this will be the combined deletion depth map, merging the for and rev observations
	for key, value := range p.DelCountF {
		dTot += value
		deletions[key] = value
	}
	for key, value := range p.DelCountR {
		dTot += value
		if _, ok := deletions[key]; ok { // if insertion on reverse strand was seen already on positive strand
			deletions[key] = deletions[key] + value //add reverse read counts to forward read counts
		} else {
			deletions[key] = value
		}
	}
	var DaKey, DaValue int
	for key, value := range deletions {
		if value > DaValue {
			DaKey = key
			DaValue = value
		}
	}
	B = N - iTot
	var noDeletionPosterior, deletionPosterior float64
	if DaValue > 0 { //if DaValue == 0, keep initial answer deletion value of 0
		noDeletionPosterior = logspace.Multiply(homozygousIndelLikelihoodExpression(B, DaValue, epsilon, homoIndelCache), haploidIndelPriorCache[0])
		deletionPosterior = logspace.Multiply(homozygousIndelLikelihoodExpression(DaValue, B, epsilon, homoIndelCache), haploidIndelPriorCache[1])
		if deletionPosterior > noDeletionPosterior {
			answer.Deletion = DaKey
		}
	}

	return answer
}

// makeHaploidBasePriorCAche is a helper function used in samAssembler before running HaploidCallFromPile.
// Constructs a matrix ([][]float64) of form answer[refBase][outputBase], where the index of both bases 0=a, 1=c, 2=g, 3=t.
// Values are log transformed.
func makeHaploidBasePriorCache(delta float64, gamma float64) [][]float64 {
	Tv := math.Log(delta / (2.0 + gamma))         //the probability of transversions
	Tr := math.Log(gamma * delta / (2.0 + gamma)) //the probability of transitions
	oneMinusDelta := math.Log(1.0 - delta)

	return [][]float64{
		{oneMinusDelta, Tv, Tr, Tv},
		{Tv, oneMinusDelta, Tv, Tr},
		{Tr, Tv, oneMinusDelta, Tv},
		{Tv, Tr, Tv, oneMinusDelta}}
}

// makeHaploidIndelPriorCache is a helper function used in before calling HaploidCallFromPile.
// Stores two values: answer[0] = p(Base), answer[1] = p(Indel).
func makeHaploidIndelPriorCache(delta float64, kappa float64) []float64 {
	return []float64{math.Log(1.0 - delta*kappa), math.Log(delta * kappa)}
}
