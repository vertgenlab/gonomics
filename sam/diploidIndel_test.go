package sam

import "testing"

var DiploidIndelCallFromPileTests = []struct {
	P                 Pile
	Delta             float64
	Epsilon           float64
	Kappa             float64
	ExpectedInsertion DiploidInsertion
	ExpectedDeletion  DiploidDeletion
}{
	{P: Pile{CountF: [13]int{30, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, //high coverage, clear complex het insertion
		InsCountF: map[string]int{"AAT": 7, "AT": 6},
		InsCountR: map[string]int{"AAT": 6, "AT": 5}},
		Delta:             0.01,
		Epsilon:           0.01,
		Kappa:             0.05,
		ExpectedInsertion: DiploidInsertion{Type: IaIb, Ia: "AAT", Ib: "AT"},
		ExpectedDeletion:  DiploidDeletion{Type: BBNoDel, Da: 0, Db: 0},
	},
	{P: Pile{CountF: [13]int{30, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, //high coverage, het simple insertion (IaB)
		InsCountF: map[string]int{"AAT": 7, "AT": 1},
		InsCountR: map[string]int{"AAT": 6}},
		Delta:             0.01,
		Epsilon:           0.01,
		Kappa:             0.05,
		ExpectedInsertion: DiploidInsertion{Type: IaB, Ia: "AAT", Ib: "AT"},
		ExpectedDeletion:  DiploidDeletion{Type: BBNoDel, Da: 0, Db: 0},
	},
	{P: Pile{CountF: [13]int{30, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, //no insertion, default parameters reject with only 4 observations (BB)
		InsCountF: map[string]int{"AAT": 1},
		InsCountR: map[string]int{"AAT": 3}},
		Delta:             0.01,
		Epsilon:           0.01,
		Kappa:             0.05,
		ExpectedInsertion: DiploidInsertion{Type: BBnoIns, Ia: "AAT", Ib: ""},
		ExpectedDeletion:  DiploidDeletion{Type: BBNoDel, Da: 0, Db: 0},
	},
	{P: Pile{CountF: [13]int{30, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, //default parameters believe with at least 5 observations (IaB)
		InsCountF: map[string]int{"AAT": 2},
		InsCountR: map[string]int{"AAT": 3}},
		Delta:             0.01,
		Epsilon:           0.01,
		Kappa:             0.05,
		ExpectedInsertion: DiploidInsertion{Type: IaB, Ia: "AAT", Ib: ""},
		ExpectedDeletion:  DiploidDeletion{Type: BBNoDel, Da: 0, Db: 0},
	},
	{P: Pile{CountF: [13]int{30, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, //magic number for a homozygous call at 30x coverage is 27. (IaIa)
		InsCountF: map[string]int{"AAT": 12},
		InsCountR: map[string]int{"AAT": 17}},
		Delta:             0.01,
		Epsilon:           0.01,
		Kappa:             0.05,
		ExpectedInsertion: DiploidInsertion{Type: IaIa, Ia: "AAT", Ib: ""},
		ExpectedDeletion:  DiploidDeletion{Type: BBNoDel, Da: 0, Db: 0},
	},
	{P: Pile{CountF: [13]int{60, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, //Homozygous, higher coverage
		InsCountF: map[string]int{"AAT": 23},
		InsCountR: map[string]int{"AAT": 34}},
		Delta:             0.01,
		Epsilon:           0.01,
		Kappa:             0.05,
		ExpectedInsertion: DiploidInsertion{Type: IaIa, Ia: "AAT", Ib: ""},
		ExpectedDeletion:  DiploidDeletion{Type: BBNoDel, Da: 0, Db: 0},
	},
	{P: Pile{CountF: [13]int{30, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, //high coverage, clear complex het deletion
		DelCountF: map[int]int{3: 7, 2: 6},
		DelCountR: map[int]int{3: 6, 2: 5}},
		Delta:             0.01,
		Epsilon:           0.01,
		Kappa:             0.05,
		ExpectedInsertion: DiploidInsertion{Type: BBnoIns, Ia: "", Ib: ""},
		ExpectedDeletion:  DiploidDeletion{Type: DaDb, Da: 3, Db: 2},
	},
	{P: Pile{CountF: [13]int{30, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, //high coverage, het simple deletion (DaB)
		DelCountF: map[int]int{3: 7, 2: 1},
		DelCountR: map[int]int{3: 6}},
		Delta:             0.01,
		Epsilon:           0.01,
		Kappa:             0.05,
		ExpectedInsertion: DiploidInsertion{Type: BBnoIns, Ia: "", Ib: ""},
		ExpectedDeletion:  DiploidDeletion{Type: DaB, Da: 3, Db: 2},
	},
}

func TestDiploidIndelCallFromPile(t *testing.T) {
	var priorCache []float64
	var emptyCache = make([][]float64, 0)
	var actualInsertion DiploidInsertion
	var actualDeletion DiploidDeletion
	for _, v := range DiploidIndelCallFromPileTests {
		priorCache = MakeDiploidIndelPriorCache(v.Kappa, v.Delta)
		actualInsertion = DiploidInsertionCallFromPile(v.P, priorCache, emptyCache, emptyCache, v.Epsilon)
		if actualInsertion != v.ExpectedInsertion {
			t.Errorf("Error in DiploidInsertionCallFromPile. Expected: %s. Found: %s.\n", diploidInsertionString(v.ExpectedInsertion), diploidInsertionString(actualInsertion))
		}
		actualDeletion = DiploidDeletionCallFromPile(v.P, priorCache, emptyCache, emptyCache, v.Epsilon)
		if actualDeletion != v.ExpectedDeletion {
			t.Errorf("Error in DiploidDeletionCallFromPile. Expected: %s. Found: %s.\n", diploidDeletionString(v.ExpectedDeletion), diploidDeletionString(actualDeletion))
		}
	}
}

var IndelLikelihoodExpressionTests = []struct {
	CorrectCount   int
	IncorrectCount int
	Epsilon        float64
	ExpectedHomo   float64
	ExpectedHetero float64
}{
	{CorrectCount: 30,
		IncorrectCount: 0,
		Epsilon:        0.01,
		ExpectedHomo:   -0.3015100756050435,
		ExpectedHetero: -20.944791671504685},
	{CorrectCount: 25,
		IncorrectCount: 4,
		Epsilon:        0.01,
		ExpectedHomo:   -21.444527862529682,
		ExpectedHetero: -38.64726252577938},
}

func TestIndelLikelihoodExpressions(t *testing.T) {
	var cache = make([][]float64, 0) //this will always be of dimension 0x0 in testing so we calculate by hand
	var actual float64
	for _, v := range IndelLikelihoodExpressionTests {
		actual = homozygousIndelLikelihoodExpression(v.CorrectCount, v.IncorrectCount, v.Epsilon, cache)
		if actual != v.ExpectedHomo {
			t.Errorf("Error in homozygousIndelLikelihoodExpression. Expected: %v. Found: %v.\n", v.ExpectedHomo, actual)
		}
		actual = heterozygousIndelLikelihoodExpression(v.CorrectCount, v.IncorrectCount, v.Epsilon, cache)
		if actual != v.ExpectedHetero {
			t.Errorf("Error in heteroIndelLikelihoodExpression. Expected: %v. Found: %v.\n", v.ExpectedHetero, actual)
		}
	}
}

var MakeDiploidIndelPriorCacheTests = []struct {
	Delta    float64
	Kappa    float64
	Expected []float64
}{
	{Kappa: 0.05,
		Delta:    0.01,
		Expected: []float64{-15.201804919084164, -14.508657738524219, -5.600902459542082, -0.0020027541739614635}},
}

func TestMakeDiploidIndelPriorCache(t *testing.T) {
	var current []float64
	var i int
	for _, v := range MakeDiploidIndelPriorCacheTests {
		current = MakeDiploidIndelPriorCache(v.Kappa, v.Delta)
		if len(current) != len(v.Expected) {
			t.Errorf("Error in makeDiploidIndelPriorCache. Answer was the wrong length.")
		}
		for i = range current {
			if current[i] != v.Expected[i] {
				t.Errorf("Error in makeDiploidIndelPriorCache. Expected: %v. Found:%v.\n", v.Expected, current)
			}
		}
	}
}
