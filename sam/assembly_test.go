package sam

import (
	"fmt"
	"github.com/vertgenlab/gonomics/dna"
	"testing"
)

var DiploidBaseCallFromPileTests = []struct {
	P        Pile
	RefBase  dna.Base
	Delta    float64
	Gamma    float64
	Epsilon  float64
	Expected DiploidBase
}{
	{P: Pile{CountF: [13]int{16, 14, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0}}, //good coverage heterozygous double mutation
		RefBase:  dna.G,
		Delta:    0.01,
		Gamma:    3,
		Epsilon:  0.01,
		Expected: AC},
	{P: Pile{CountF: [13]int{4, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0}}, //little data, only confident enough to call one mutation
		RefBase:  dna.G,
		Delta:    0.01,
		Gamma:    3,
		Epsilon:  0.01,
		Expected: AG},
	{P: Pile{CountF: [13]int{4, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0}}, //changing expected divergence rate changes pile outcome
		RefBase:  dna.G,
		Delta:    0.5,
		Gamma:    3,
		Epsilon:  0.01,
		Expected: AT},
	{P: Pile{CountF: [13]int{4, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0}}, //changing expected error rate changes pile outcome.
		RefBase:  dna.G,
		Delta:    0.01,
		Gamma:    3,
		Epsilon:  0.0001,
		Expected: AT},
	{P: Pile{CountF: [13]int{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}}, //test empty pile, returns homozygous reference
		RefBase:  dna.G,
		Delta:    0.01,
		Gamma:    3,
		Epsilon:  0.01,
		Expected: GG},
	{P: Pile{CountF: [13]int{16, 450, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0}}, //stress test, huge pile
		RefBase:  dna.G,
		Delta:    0.01,
		Gamma:    3,
		Epsilon:  0.01,
		Expected: CC},
	{P: Pile{CountF: [13]int{16, 14, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0}}, //reference N returns NN.
		RefBase:  dna.N,
		Delta:    0.01,
		Gamma:    3,
		Epsilon:  0.01,
		Expected: NN},
	{P: Pile{CountF: [13]int{16, 1, 0, 7, 0, 0, 0, 0, 0, 0, 0, 0, 0}}, //messy pile, high error rate
		RefBase:  dna.C,
		Delta:    0.1,
		Gamma:    3,
		Epsilon:  0.01,
		Expected: AT},
}

func TestDiploidBaseCallFromPile(t *testing.T) {
	var actual DiploidBase
	var priorCache [][]float64
	var heterozycousCache = make([][]float64, 0)
	var homozygousCache = make([][]float64, 0)

	for _, v := range DiploidBaseCallFromPileTests {
		priorCache = makePriorCache(v.Delta, v.Gamma)
		actual = DiploidBaseCallFromPile(v.P, v.RefBase, priorCache, homozygousCache, heterozycousCache, v.Epsilon)
		if actual != v.Expected {
			t.Errorf("Error in DiploidBaseCallFromPile. Expected: %s. Observed: %s.", diploidBaseString(v.Expected), diploidBaseString(actual))
		}
	}
}

var LikelihoodExpressionTests = []struct {
	CorrectCount   int
	IncorrectCount int
	Epsilon        float64
	Cache          [][]float64 //this will always be o dimension 0x0 in testing so we calculate by hand
	ExpectedHomo   float64
	ExpectedHetero float64
}{
	{
		CorrectCount:   26,
		IncorrectCount: 3,
		Epsilon:        0.01,
		Cache:          [][]float64{},
		ExpectedHomo:   -17.37265615615964,
		ExpectedHetero: -35.3070878104479,
	},
	{
		CorrectCount:   14,
		IncorrectCount: 16,
		Epsilon:        0.01,
		Cache:          [][]float64{},
		ExpectedHomo:   -91.40122429644823,
		ExpectedHetero: -101.0582259564496,
	},
}

func TestLikelihoodExpressions(t *testing.T) {
	var actual float64
	for _, v := range LikelihoodExpressionTests {
		actual = homozygousLikelihoodExpression(v.CorrectCount, v.IncorrectCount, v.Epsilon, v.Cache)
		if actual != v.ExpectedHomo {
			t.Errorf("Error in homozygousLikelihoodExpression. Expected: %v. Found:%v.", v.ExpectedHomo, actual)
		}
		actual = heterozygousLikelihoodExpression(v.CorrectCount, v.IncorrectCount, v.Epsilon, v.Cache)
		if actual != v.ExpectedHetero {
			t.Errorf("Error in heterozygousLikelihoodExpression. Expected: %v. Found: %v.", v.ExpectedHetero, actual)
		}
	}
}

var MakePileDiploidPriorProbabilityCacheTests = []struct {
	BranchLength   float64
	TransitionBias float64
	Expected       [][]float64
}{
	{BranchLength: 0.01,
		TransitionBias: 3,
		Expected: [][]float64{[]float64{-0.02010067170700291, -5.531511253715748, -4.432898965047638, -5.531511253715748, -12.429216196844383, -11.330603908176274, -11.736069016284437, -10.231991619508165, -11.330603908176274, -12.429216196844383},
			[]float64{-12.429216196844383, -5.531511253715748, -11.736069016284437, -11.330603908176274, -0.02010067170700291, -5.531511253715748, -4.432898965047638, -12.429216196844383, -11.330603908176274, -10.231991619508165},
			[]float64{-10.231991619508165, -11.330603908176274, -4.432898965047638, -11.330603908176274, -12.429216196844383, -5.531511253715748, -11.736069016284437, -0.02010067170700291, -5.531511253715748, -12.429216196844383},
			[]float64{-12.429216196844383, -11.330603908176274, -11.736069016284437, -5.531511253715748, -10.231991619508165, -11.330603908176274, -4.432898965047638, -11.736069016284437, -5.531511253715748, -0.02010067170700291}},
	},
}

func TestMakePileDiploidPriorProbabilityCache(t *testing.T) {
	var current [][]float64
	for _, v := range MakePileDiploidPriorProbabilityCacheTests {
		current = makePriorCache(v.BranchLength, v.TransitionBias)
		if !equalMatrix(current, v.Expected) {
			fmt.Println(current)
			t.Errorf("Error in generating pileup DiploidBase prior probability cache.")
		}
	}
}

func equalMatrix(a [][]float64, b [][]float64) bool {
	var i, j int
	if len(a) != len(b) {
		return false
	}
	for i = range a {
		if len(a[i]) != len(b[i]) {
			return false
		}
		for j = range a[i] {
			if a[i][j] != b[i][j] {
				return false
			}
		}
	}
	return true
}
