package sam

import (
	"testing"

	"github.com/vertgenlab/gonomics/dna"
)

var HaploidCallFromPileTests = []struct {
	P        Pile
	RefBase  dna.Base
	Gamma    float64
	Delta    float64
	Epsilon  float64
	Kappa    float64
	Lambda   float64
	Expected HaploidCall
}{
	{P: Pile{CountF: [13]int{30, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, //messy start insertion
		InsCountF: map[string]int{"AAT": 7, "AT": 6},
		InsCountR: map[string]int{"AAT": 6, "AT": 5}},
		RefBase:  dna.C,
		Gamma:    3,
		Delta:    0.01,
		Epsilon:  0.01,
		Kappa:    0.05,
		Expected: HaploidCall{Base: dna.A, Insertion: "AAT", Deletion: 0},
	},
	{P: Pile{CountF: [13]int{30, 0, 13, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, //messy start insertion
		InsCountF: map[string]int{"AAT": 2, "AT": 6},
		InsCountR: map[string]int{"AAT": 3, "AT": 5},
		DelCountF: map[int]int{2: 3, 5: 29}},
		RefBase:  dna.A,
		Gamma:    3,
		Delta:    0.01,
		Epsilon:  0.01,
		Kappa:    0.05,
		Expected: HaploidCall{Base: dna.A, Insertion: "", Deletion: 5},
	},
	{P: Pile{CountF: [13]int{30, 0, 5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, //messy start insertion
		InsCountF: map[string]int{"AAT": 2, "AT": 6},
		InsCountR: map[string]int{},
		DelCountF: map[int]int{2: 30}},
		RefBase:  dna.A,
		Gamma:    3,
		Delta:    0.1,
		Epsilon:  0.01,
		Kappa:    0.5,
		Expected: HaploidCall{Base: dna.A, Insertion: "", Deletion: 2},
	},
	{P: Pile{CountF: [13]int{30, 13, 5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, //high lambda, G -> A common
		InsCountF: map[string]int{"AAT": 2, "AT": 6},
		InsCountR: map[string]int{},
		DelCountF: map[int]int{2: 46}},
		RefBase:  dna.A,
		Gamma:    3,
		Delta:    0.1,
		Epsilon:  0.01,
		Kappa:    0.5,
		Lambda:   0.5,
		Expected: HaploidCall{Base: dna.G, Insertion: "", Deletion: 2},
	},
}

func TestHaploidCallFromPile(t *testing.T) {
	var homoBaseCache = make([][]float64, 0)
	var heteroBaseCache = make([][]float64, 0)
	var homoIndelCache = make([][]float64, 0)
	var haploidBasePriorCache [][]float64
	var haploidIndelPriorCache []float64
	var observed HaploidCall
	var ancientCache = AncientLikelihoodCache{
		EpsilonOverThreePlusLambdaOverTwo:                []float64{},
		OneMinusEpsilon:                                  []float64{},
		OneMinusEpsilonMinusLambda:                       []float64{},
		EpsilonOverThreePlusLambda:                       []float64{},
		PointFiveMinusEpsilonOverThreeMinusLambdaOverTwo: []float64{},
		EpsilonOverThree:                                 []float64{},
		PointFiveMinusEpsilonOverThreePlusLambdaOverTwo:  []float64{},
		PointFiveMinusEpsilonOverThree:                   []float64{},
	}
	for _, v := range HaploidCallFromPileTests {
		haploidBasePriorCache = MakeHaploidBasePriorCache(v.Delta, v.Gamma)
		haploidIndelPriorCache = MakeHaploidIndelPriorCache(v.Delta, v.Kappa)
		observed = HaploidCallFromPile(v.P, v.RefBase, v.Epsilon, v.Lambda, haploidBasePriorCache, haploidIndelPriorCache, homoBaseCache, heteroBaseCache, homoIndelCache, ancientCache)
		if observed.Base != v.Expected.Base {
			t.Errorf("Error in HaploidCallFromPile. Observed base: %v not equal to expected base: %v.\n", dna.BaseToString(observed.Base), dna.BaseToString(v.Expected.Base))
		}
		if observed.Insertion != v.Expected.Insertion {
			t.Errorf("Error in HaploidCallFromPile. Observed insertion: %v not equal to expected insertion: %v.\n", observed.Insertion, v.Expected.Insertion)
		}
		if observed.Deletion != v.Expected.Deletion {
			t.Errorf("Error in HaploidCallFromPile. Observed deletion: %v not equal to expected deletion: %v.\n", observed.Deletion, v.Expected.Deletion)
		}
	}
}
