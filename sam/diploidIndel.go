package sam

import (
	"fmt"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/numbers"
	"github.com/vertgenlab/gonomics/numbers/logspace"
	"log"
	"math"
)

// InsertionType encodes the insertion genotype state, which can be one of the four constant values explained below.
type InsertionType byte

const (
	IaIa    InsertionType = 0 // for insertionA insertionA, a homozygous insertion
	IaIb    InsertionType = 1 // for insertionA insertionB, or complex insertion
	IaB     InsertionType = 2 // for insertion base, or heterozygous insertion
	BBnoIns InsertionType = 3 // for base base, or no insertion
)

// DiploidInsertion is a minimal internal format to encode Insertion variant genotypes
type DiploidInsertion struct {
	Type InsertionType
	Ia   string
	Ib   string
}

// diploidInsertionString formats a DiploidInsertion struct as a string for debugging.
func diploidInsertionString(i DiploidInsertion) string {
	switch i.Type {
	case IaIa:
		return fmt.Sprintf("IaIa. Ia: %s.\n", i.Ia)
	case IaIb:
		return fmt.Sprintf("IaIb. Ia: %s. Ib: %s.\n", i.Ia, i.Ib)
	case IaB:
		return fmt.Sprintf("IaB. Ia: %s.\n", i.Ia)
	case BBnoIns:
		return fmt.Sprintf("BB.\n")
	}
	log.Fatalf("InsertionType not recognized.")
	return ""
}

// DiploidInsertionCallFromPile returns a DiploidInsertion struct from an input Pile after making an insertion variant call.
// Takes in a cache for the prior distribution values, and two likelihood caches.
// Epsilon defines the misclassification rate parameter.
func DiploidInsertionCallFromPile(p Pile, priorCache []float64, homozygousIndelCache [][]float64, heterozygousIndelCache [][]float64, epsilon float64) DiploidInsertion {
	var iTot int

	var aCount = p.CountF[dna.A] + p.CountR[dna.A]
	var cCount = p.CountF[dna.C] + p.CountR[dna.C]
	var gCount = p.CountF[dna.G] + p.CountR[dna.G]
	var tCount = p.CountF[dna.T] + p.CountR[dna.T]
	var N = aCount + cCount + gCount + tCount //this is the base read depth at the pile.

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

	var IaKey, IbKey string
	var IaValue, IbValue int

	for key, value := range insertions {
		if value > IaValue {
			IbKey = IaKey
			IbValue = IaValue
			IaKey = key
			IaValue = value
		} else if value > IbValue {
			IbKey = key
			IbValue = value
		}
	}

	if IaValue < 1 { //if no insertion counts were observed, return BaseBase
		return DiploidInsertion{Type: BBnoIns, Ia: "", Ib: ""}
	}

	//DEBUG: fmt.Printf("N: %v. dTot: %v. DaValue: %v. DbValue: %v.\n", N, iTot, IaValue, IbValue)

	var B int = N - iTot
	var answer []DiploidInsertion = []DiploidInsertion{DiploidInsertion{Type: BBnoIns, Ia: IaKey, Ib: IbKey}}
	answerPosterior := logspace.Multiply(homozygousIndelLikelihoodExpression(B, IaValue+IbValue, epsilon, homozygousIndelCache), priorCache[BBnoIns])
	//DEBUG: fmt.Printf("BB Insertion Posterior: %v.\n", answerPosterior)

	var currentPosterior float64
	//now we check the other genotypes
	//IaIa
	currentPosterior = logspace.Multiply(homozygousIndelLikelihoodExpression(IaValue, B+IbValue, epsilon, homozygousIndelCache), priorCache[IaIa])
	//DEBUG: fmt.Printf("IaIa Insertion Posterior: %v.\n", currentPosterior)
	if currentPosterior > answerPosterior {
		answer = answer[:1] //clear ties
		answer[0].Type = IaIa
		answerPosterior = currentPosterior
	} else if currentPosterior == answerPosterior {
		answer = append(answer, DiploidInsertion{Type: IaIa, Ia: IaKey, Ib: IbKey})
	}
	//IaIb
	currentPosterior = logspace.Multiply(heterozygousIndelLikelihoodExpression(IaValue+IbValue, B, epsilon, heterozygousIndelCache), priorCache[IaIb])
	//DEBUG: fmt.Printf("IaIb Insertion Posterior: %v.\n", currentPosterior)
	if currentPosterior > answerPosterior {
		answer = answer[:1] //clear ties
		answer[0].Type = IaIb
		answerPosterior = currentPosterior
	} else if currentPosterior == answerPosterior {
		answer = append(answer, DiploidInsertion{Type: IaIb, Ia: IaKey, Ib: IbKey})
	}
	//IaBase
	currentPosterior = logspace.Multiply(heterozygousIndelLikelihoodExpression(IaValue+B, IbValue, epsilon, heterozygousIndelCache), priorCache[IaB])
	//DEBUG: fmt.Printf("IaB Insertion Posterior: %v.\n", currentPosterior)
	if currentPosterior > answerPosterior {
		answer = answer[:1] //clear ties
		answer[0].Type = IaB
		answerPosterior = currentPosterior
	} else if currentPosterior == answerPosterior {
		answer = append(answer, DiploidInsertion{Type: IaB, Ia: IaKey, Ib: IbKey})
	}

	return answer[numbers.RandIntInRange(0, len(answer))] //if two insertion genotypes have the same posterior probability, pick one at random
}

// DeletionType encodes the deletion genotype state, which can be one of the four constant values explained below.
type DeletionType byte

const (
	DaDa    DeletionType = 0
	DaDb    DeletionType = 1
	DaB     DeletionType = 2
	BBNoDel DeletionType = 3
)

// DiploidDeletion is a minimal internal format to encode Deletion variant genotypes
type DiploidDeletion struct {
	Type DeletionType
	Da   int
	Db   int
}

// diploidDeletionString formats a DiploidDeletion struct as a string for debugging.
func diploidDeletionString(i DiploidDeletion) string {
	switch i.Type {
	case DaDa:
		return fmt.Sprintf("DaDa. Da: %v.\n", i.Da)
	case DaDb:
		return fmt.Sprintf("DaDb. Da: %v. Db: %v.\n", i.Da, i.Db)
	case DaB:
		return fmt.Sprintf("DaB. Da: %v.\n", i.Da)
	case BBNoDel:
		return fmt.Sprintf("BB.\n")
	}
	log.Fatalf("DeletionType not recognized.")
	return ""
}

// DiploidDeletionCallFromPile returns a DiploidDeletion struct from an input Pile after making a deletion variant call.
// Takes in a cache for the prior distribution values, and two likelihood caches.
// Epsilon defines the misclassification rate parameter.
func DiploidDeletionCallFromPile(p Pile, priorCache []float64, homozygousIndelCache [][]float64, heterozygousIndelCache [][]float64, epsilon float64) DiploidDeletion {
	var dTot int

	var aCount = p.CountF[dna.A] + p.CountR[dna.A]
	var cCount = p.CountF[dna.C] + p.CountR[dna.C]
	var gCount = p.CountF[dna.G] + p.CountR[dna.G]
	var tCount = p.CountF[dna.T] + p.CountR[dna.T]
	var N = aCount + cCount + gCount + tCount //this is the base read depth at the pile.

	var deletions = make(map[int]int, 0) //this will be the combined deletion depth map, merging the for and rev observations
	for key, value := range p.DelCountF {
		dTot += value
		deletions[key] = value
	}
	for key, value := range p.DelCountR {
		dTot += value
		if _, ok := deletions[key]; ok { // if deletion on reverse strand was seen already on positive strand
			deletions[key] = deletions[key] + value //add reverse read counts to forward read counts
		} else {
			deletions[key] = value
		}
	}

	var DaKey, DbKey, DaValue, DbValue int

	for key, value := range deletions {
		if value > DaValue {
			DbKey = DaKey
			DbValue = DaValue
			DaKey = key
			DaValue = value
		} else if value > DbValue {
			DbKey = key
			DbValue = value
		}
	}

	if DaValue < 1 { //if no deletion counts were observed, return BaseBase
		return DiploidDeletion{Type: BBNoDel, Da: 0, Db: 0}
	}

	//DEBUG: fmt.Printf("N: %v. dTot: %v. DaValue: %v. DbValue: %v.\n", N, dTot, DaValue, DbValue)

	var B int = N - dTot

	//set default return to no deletion
	var answer []DiploidDeletion = []DiploidDeletion{DiploidDeletion{Type: BBNoDel, Da: DaKey, Db: DbKey}}
	answerPosterior := logspace.Multiply(homozygousIndelLikelihoodExpression(B, DaValue+DbValue, epsilon, homozygousIndelCache), priorCache[BBNoDel])

	var currentPosterior float64
	//now we check the other genotypes
	//DaDa
	currentPosterior = logspace.Multiply(homozygousIndelLikelihoodExpression(DaValue, B+DbValue, epsilon, homozygousIndelCache), priorCache[DaDa])
	if currentPosterior > answerPosterior {
		answer = answer[:1] //clear ties
		answer[0].Type = DaDa
		answerPosterior = currentPosterior
	} else if currentPosterior == answerPosterior {
		answer = append(answer, DiploidDeletion{Type: DaDa, Da: DaKey, Db: DbKey})
	}
	//DaDb
	currentPosterior = logspace.Multiply(heterozygousIndelLikelihoodExpression(DaValue+DbValue, B, epsilon, heterozygousIndelCache), priorCache[DaDb])
	if currentPosterior > answerPosterior {
		answer = answer[:1] //clear ties
		answer[0].Type = DaDb
		answerPosterior = currentPosterior
	} else if currentPosterior == answerPosterior {
		answer = append(answer, DiploidDeletion{Type: DaDb, Da: DaKey, Db: DbKey})
	}
	//IaBase
	currentPosterior = logspace.Multiply(heterozygousIndelLikelihoodExpression(DaValue+B, DbValue, epsilon, heterozygousIndelCache), priorCache[DaB])
	if currentPosterior > answerPosterior {
		answer = answer[:1] //clear ties
		answer[0].Type = DaB
		answerPosterior = currentPosterior
	} else if currentPosterior == answerPosterior {
		answer = append(answer, DiploidDeletion{Type: DaB, Da: DaKey, Db: DbKey})
	}

	return answer[numbers.RandIntInRange(0, len(answer))] //if two insertion genotypes have the same posterior probability, pick one at random
}

// homozygousIndelLikelihoodExpression is a helper function of DiploidInsertionCallFromPile and DiploidDeletionCallFromPile
// and calculates the multinomial expression for homozygous INDEL genotypes.
func homozygousIndelLikelihoodExpression(correctCount int, incorrectCount int, epsilon float64, homozygousIndelCache [][]float64) float64 {
	//DEBUG: fmt.Printf("Correct: %v. Incorrect: %v.\n", correctCount, incorrectCount)
	if correctCount < len(homozygousIndelCache) && incorrectCount < len(homozygousIndelCache[correctCount]) { //if the indel coverage is within the cache bounds
		if homozygousIndelCache[correctCount][incorrectCount] != 0 {
			return homozygousIndelCache[correctCount][incorrectCount]
		} else {
			s := logspace.Pow(math.Log(1.0-epsilon), float64(correctCount))
			f := logspace.Pow(math.Log(epsilon/2.0), float64(incorrectCount))
			homozygousIndelCache[correctCount][incorrectCount] = logspace.Multiply(s, f)
			return homozygousIndelCache[correctCount][incorrectCount]
		}
	} else {
		s := logspace.Pow(math.Log(1.0-epsilon), float64(correctCount))
		f := logspace.Pow(math.Log(epsilon/2.0), float64(incorrectCount))
		return logspace.Multiply(s, f)
	}
}

// heterozygousIndelLikelihoodExpression is a helper function of DiploidInsertionCallFromPile and DiploidDeletionCallFromPile
// and calculates the multinomial expression for heterozygous INDEL genotypes.
func heterozygousIndelLikelihoodExpression(correctCount int, incorrectCount int, epsilon float64, heterozygousIndelCache [][]float64) float64 {
	//DEBUG: fmt.Printf("Correct: %v. Incorrect: %v.\n", correctCount, incorrectCount)
	if correctCount < len(heterozygousIndelCache) && incorrectCount < len(heterozygousIndelCache[correctCount]) { //if the indel coverage is within the cache bounds
		if heterozygousIndelCache[correctCount][incorrectCount] != 0 {
			return heterozygousIndelCache[correctCount][incorrectCount]
		} else {
			s := logspace.Pow(math.Log(0.5-(epsilon/4.0)), float64(correctCount))
			f := logspace.Pow(math.Log(epsilon/2.0), float64(incorrectCount))
			heterozygousIndelCache[correctCount][incorrectCount] = logspace.Multiply(s, f)
			return heterozygousIndelCache[correctCount][incorrectCount]
		}
	} else {
		s := logspace.Pow(math.Log(0.5-(epsilon/4.0)), float64(correctCount))
		f := logspace.Pow(math.Log(epsilon/2.0), float64(incorrectCount))
		return logspace.Multiply(s, f)
	}
}

// makeDiploidIndelPriorCache is a helper function used in samAssembler before running
// DiploidInsertionCallFromPile and DiploidDeletionCallFromPile. Constructs a []float64, where the index corresponds to InsertionType / DeletionType,
// and the value refers to the log transformed prior probability density.
// parameterized on delta, the expected divergence rate, and kappa, the proportion of mutations expected to be INDELs.
func makeDiploidIndelPriorCache(kappa float64, delta float64) []float64 {
	kd := logspace.Multiply(math.Log(kappa), math.Log(delta))
	kdSquared := logspace.Pow(kd, 2)
	pBaseBase := math.Log(1 - 4*kappa*delta - 3*(kappa*kappa*delta*delta))
	return []float64{kdSquared, logspace.Multiply(math.Log(2), kdSquared), logspace.Multiply(2, kd), pBaseBase}
}
