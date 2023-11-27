package sam

import (
	"testing"

	"github.com/vertgenlab/gonomics/dna"
)

var DiploidBaseCallFromPileTests = []struct {
	P        Pile
	RefBase  dna.Base
	Delta    float64
	Gamma    float64
	Epsilon  float64
	Lambda   float64
	Expected DiploidBase
}{
	{P: Pile{CountF: [13]int{16, 14, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0}}, //good coverage heterozygous double mutation
		RefBase:  dna.G,
		Delta:    0.01,
		Gamma:    3,
		Epsilon:  0.01,
		Lambda:   0,
		Expected: AC},
	{P: Pile{CountF: [13]int{4, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0}}, //little data, only confident enough to call one mutation
		RefBase:  dna.G,
		Delta:    0.01,
		Gamma:    3,
		Epsilon:  0.01,
		Lambda:   0,
		Expected: AG},
	{P: Pile{CountF: [13]int{4, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0}}, //changing expected divergence rate changes pile outcome
		RefBase:  dna.G,
		Delta:    0.5,
		Gamma:    3,
		Epsilon:  0.01,
		Lambda:   0,
		Expected: AT},
	{P: Pile{CountF: [13]int{4, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0}}, //changing expected error rate changes pile outcome.
		RefBase:  dna.G,
		Delta:    0.01,
		Gamma:    3,
		Epsilon:  0.0001,
		Lambda:   0,
		Expected: AT},
	{P: Pile{CountF: [13]int{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}}, //test empty pile, returns homozygous reference
		RefBase:  dna.G,
		Delta:    0.01,
		Gamma:    3,
		Epsilon:  0.01,
		Lambda:   0,
		Expected: GG},
	{P: Pile{CountF: [13]int{16, 450, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0}}, //stress test, huge pile
		RefBase:  dna.G,
		Delta:    0.01,
		Gamma:    3,
		Epsilon:  0.01,
		Lambda:   0,
		Expected: CC},
	{P: Pile{CountF: [13]int{16, 14, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0}}, //reference N returns NN.
		RefBase:  dna.N,
		Delta:    0.01,
		Gamma:    3,
		Epsilon:  0.01,
		Lambda:   0,
		Expected: NN},
	{P: Pile{CountF: [13]int{16, 1, 0, 7, 0, 0, 0, 0, 0, 0, 0, 0, 0}}, //messy pile, high error rate
		RefBase:  dna.C,
		Delta:    0.1,
		Gamma:    3,
		Epsilon:  0.01,
		Lambda:   0,
		Expected: AT},
	{P: Pile{CountF: [13]int{0, 61, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}}, //homo transversion
		RefBase:  dna.T,
		Delta:    0.1,
		Gamma:    3,
		Epsilon:  0.01,
		Lambda:   0,
		Expected: CC},
	{P: Pile{CountF: [13]int{16, 14, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0}}, //first lambda test
		RefBase:  dna.G,
		Delta:    0.01,
		Gamma:    3,
		Epsilon:  0.01,
		Lambda:   0.05,
		Expected: AC},
	{P: Pile{CountF: [13]int{16, 4, 0, 7, 0, 0, 0, 0, 0, 0, 0, 0, 0}}, //messy pile, high lambda
		RefBase:  dna.C,
		Delta:    0.1,
		Gamma:    3,
		Epsilon:  0.01,
		Lambda:   0.2,
		Expected: AC},
}

func TestDiploidBaseCallFromPile(t *testing.T) {
	var actual DiploidBase
	var priorCache [][]float64
	var heterozygousCache = make([][]float64, 0)
	var homozygousCache = make([][]float64, 0)
	var cache = make([]float64, 0) // this will serve as an empty cache for all expression tests
	var ancientCache = AncientLikelihoodCache{
		EpsilonOverThree:                                 cache,
		OneMinusEpsilon:                                  cache,
		OneMinusEpsilonMinusLambda:                       cache,
		EpsilonOverThreePlusLambda:                       cache,
		PointFiveMinusEpsilonOverThree:                   cache,
		EpsilonOverThreePlusLambdaOverTwo:                cache,
		PointFiveMinusEpsilonOverThreePlusLambdaOverTwo:  cache,
		PointFiveMinusEpsilonOverThreeMinusLambdaOverTwo: cache,
	}

	for _, v := range DiploidBaseCallFromPileTests {
		priorCache = MakeDiploidBasePriorCache(v.Delta, v.Gamma)
		actual = DiploidBaseCallFromPile(v.P, v.RefBase, priorCache, homozygousCache, heterozygousCache, ancientCache, v.Epsilon, v.Lambda)
		if actual != v.Expected {
			t.Errorf("Error in DiploidBaseCallFromPile. Expected: %s. Observed: %s.", DiploidBaseString(v.Expected), DiploidBaseString(actual))
		}
	}
}

var LikelihoodExpressionTests = []struct {
	CorrectCount   int
	IncorrectCount int
	Epsilon        float64
	ExpectedHomo   float64
	ExpectedHetero float64
}{
	{
		CorrectCount:   26,
		IncorrectCount: 3,
		Epsilon:        0.01,
		ExpectedHomo:   -17.37265615615964,
		ExpectedHetero: -35.3070878104479,
	},
	{
		CorrectCount:   14,
		IncorrectCount: 16,
		Epsilon:        0.01,
		ExpectedHomo:   -91.40122429644823,
		ExpectedHetero: -101.0582259564496,
	},
}

func TestBaseLikelihoodExpressions(t *testing.T) {
	var actual float64
	var cache = make([][]float64, 0) //this will always be of dimension 0x0 in testing such that we calculate by hand
	for _, v := range LikelihoodExpressionTests {
		actual = homozygousLikelihoodExpression(v.CorrectCount, v.IncorrectCount, v.Epsilon, cache)
		if actual != v.ExpectedHomo {
			t.Errorf("Error in homozygousLikelihoodExpression. Expected: %v. Found:%v.", v.ExpectedHomo, actual)
		}
		actual = heterozygousLikelihoodExpression(v.CorrectCount, v.IncorrectCount, v.Epsilon, cache)
		if actual != v.ExpectedHetero {
			t.Errorf("Error in heterozygousLikelihoodExpression. Expected: %v. Found: %v.", v.ExpectedHetero, actual)
		}
	}
}

var MakeDiploidBaseEmpiricalPriorCacheTests = []struct {
	InFile    string
	Expected  [][]float64
	Precision float64
}{
	{InFile: "testdata/samAssemblerPrior.txt",
		Expected: [][]float64{
			{0.9126446419587451, 0.017122253899044093, 0.053010229750125745, 0.013935938286097592, 0.0003521717256414554, 0.0003521717256414554, 0.00018447090390742905, 0.0013583766560456138, 0.0008552741908435346, 0.00018447090390742905},
			{0.00026842362127867255, 0.02027818448023426, 2.440214738897023e-05, 0.0005124450951683749, 0.9019277696437287, 0.01612981942410932, 0.059077598828696926, 2.440214738897023e-05, 0.00026842362127867255, 0.001488530990727184},
			{0.0005128205128205129, 0.0005128205128205129, 0.05326007326007327, 2.4420024420024423e-05, 2.4420024420024423e-05, 0.016874236874236875, 0.0002686202686202687, 0.9106471306471308, 0.017606837606837608, 0.0002686202686202687},
			{1.709986320109439e-05, 0.0005300957592339261, 0.00018809849521203833, 0.01591997264021888, 0.0005300957592339261, 0.0003590971272229822, 0.04704172366621068, 1.709986320109439e-05, 0.018313953488372094, 0.9170827633378934},
		},
		Precision: 1e-6,
	},
}

func TestMakeDiploidBaseEmpiricalPriorCache(t *testing.T) {
	var current [][]float64
	for _, v := range MakeDiploidBaseEmpiricalPriorCacheTests {
		current = MakeDiploidBaseEmpiricalPriorCache(v.InFile)
		if !equalMatrix(current, v.Expected, v.Precision) {
			t.Errorf("Error in parsing an empirical prior file for samAssembler.\n")
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
		Expected: [][]float64{{-0.02010067170700291, -5.531511253715748, -4.432898965047638, -5.531511253715748, -12.429216196844383, -11.330603908176274, -11.736069016284437, -10.231991619508165, -11.330603908176274, -12.429216196844383},
			{-12.429216196844383, -5.531511253715748, -11.736069016284437, -11.330603908176274, -0.02010067170700291, -5.531511253715748, -4.432898965047638, -12.429216196844383, -11.330603908176274, -10.231991619508165},
			{-10.231991619508165, -11.330603908176274, -4.432898965047638, -11.330603908176274, -12.429216196844383, -5.531511253715748, -11.736069016284437, -0.02010067170700291, -5.531511253715748, -12.429216196844383},
			{-12.429216196844383, -11.330603908176274, -11.736069016284437, -5.531511253715748, -10.231991619508165, -11.330603908176274, -4.432898965047638, -11.736069016284437, -5.531511253715748, -0.02010067170700291}},
	},
}

func TestMakeDiploidBasePriorCache(t *testing.T) {
	var current [][]float64
	for _, v := range MakePileDiploidPriorProbabilityCacheTests {
		current = MakeDiploidBasePriorCache(v.BranchLength, v.TransitionBias)
		if !equalMatrix(current, v.Expected, 1e-6) {
			t.Errorf("Error in generating pileup DiploidBase prior probability cache.")
		}
	}
}

func equalMatrix(a [][]float64, b [][]float64, precision float64) bool {
	var i, j int
	if len(a) != len(b) {
		return false
	}
	for i = range a {
		if len(a[i]) != len(b[i]) {
			return false
		}
		for j = range a[i] {
			if (a[i][j]-b[i][j])/a[i][j] > precision {
				return false
			}
		}
	}
	return true
}
