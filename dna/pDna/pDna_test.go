package pDna

import (
	"math/rand"
	"testing"
	// "fmt"
)

var EqualBaseTests = []struct {
	P         Float32Base
	Q         Float32Base
	Precision float32
	Expected  bool
}{
	{P: Float32Base{
		A: 0.3,
		C: 0.2,
		G: 0.5,
		T: 0,
	},
		Q: Float32Base{
			A: 0.3,
			C: 0.2,
			G: 0.5,
			T: 0,
		},
		Precision: 1e-6,
		Expected:  true,
	},
	{P: Float32Base{
		A: 0.3,
		C: 0.2,
		G: 0.5,
		T: 0,
	},
		Q: Float32Base{
			A: 0.3,
			C: 0.1,
			G: 0.5,
			T: 0.1,
		},
		Precision: 1e-6,
		Expected:  false,
	},
}

func TestEqualBase(t *testing.T) {
	var observed bool
	for _, v := range EqualBaseTests {
		observed = EqualBase(v.P, v.Q, v.Precision)
		if observed != v.Expected {
			t.Errorf("Error: in pDna, EqualBase. Output was not as expected.\n")
		}
	}
}

var EntropyTests = []struct {
	Base     Float32Base
	Expected float64
}{
	{Base: Float32Base{
		A: 1,
		C: 0,
		G: 0,
		T: 0,
	},
		Expected: 0,
	},
	{Base: Float32Base{
		A: 0.25,
		C: 0.25,
		G: 0.25,
		T: 0.25,
	},
		Expected: 2,
	},
	{Base: Float32Base{
		A: 0.5,
		C: 0.25,
		G: 0.25,
		T: 0,
	},
		Expected: 1.5,
	},
}

func TestEntropy(t *testing.T) {
	var observed float64
	for _, v := range EntropyTests {
		observed = Entropy(v.Base)
		if observed != v.Expected {
			t.Errorf("Error: pDna base entropy is not as expected. Observed: %v. Expected: %v.\n", observed, v.Expected)
		}
	}
}

var ScaleTests = []struct {
	Base       Float32Base
	Multiplier float32
	Expected   Float32Base
}{
	{Base: Float32Base{
		A: 1,
		C: 0,
		G: 0,
		T: 0,
	},
		Multiplier: .3,
		Expected: Float32Base{
			A: .3,
			C: 0,
			G: 0,
			T: 0,
		},
	},
	{Base: Float32Base{
		A: 0.25,
		C: 0.25,
		G: 0.25,
		T: 0.25,
	},
		Multiplier: .2,
		Expected: Float32Base{
			A: .05,
			C: .05,
			G: .05,
			T: .05,
		},
	},
	{Base: Float32Base{
		A: 0.5,
		C: 0.25,
		G: 0.25,
		T: 0,
	},
		Multiplier: .9,
		Expected: Float32Base{
			A: 0.45,
			C: 0.225,
			G: 0.225,
			T: 0,
		},
	},
}

func TestScale(t *testing.T) {
	var observed Float32Base
	for _, v := range ScaleTests {
		observed = Scale(v.Base, v.Multiplier)
		if observed != v.Expected {
			t.Errorf("Error: scaled pDna base is not as expected. Observed: %v. Expected: %v.\n", observed, v.Expected)
		}
	}
}

var SumTests = []struct {
	Base1    Float32Base
	Base2    Float32Base
	Expected Float32Base
}{
	{Base1: Float32Base{
		A: 1,
		C: 0,
		G: .25,
		T: 0.5,
	},
		Base2: Float32Base{
			A: 1,
			C: 0,
			G: .75,
			T: .5,
		},
		Expected: Float32Base{
			A: 2,
			C: 0,
			G: 1,
			T: 1,
		},
	},
}

func TestSum(t *testing.T) {
	var observed Float32Base
	for _, v := range SumTests {
		observed = Sum(v.Base1, v.Base2)
		if observed != v.Expected {
			t.Errorf("Error: scaled pDna base is not as expected. Observed: %v. Expected: %v.\n", observed, v.Expected)
		}
	}
}

var SumsToOneTests = []struct {
	Base      Float32Base
	Precision float32
	Expected  bool
}{
	{Base: Float32Base{
		A: 0.25,
		C: 0.25,
		G: 0.25,
		T: 0.25,
	},
		Precision: 1e-3,
		Expected:  true},
	{Base: Float32Base{
		A: 0.3345,
		C: 0.3362,
		G: 0.33,
		T: 0.000009,
	},
		Precision: 1e-4,
		Expected:  false},
}

func TestSumsToOne(t *testing.T) {
	for _, v := range SumsToOneTests {
		if SumsToOne(v.Base, v.Precision) != v.Expected {
			t.Errorf("Error: pDna sumsToOne not as expected.\n")
		}
	}
}

var RandBaseTests = []struct {
	SeedSet    bool
	SetSeed    int64
	RandSource *rand.Rand
	Expected   Float32Base
	Precision  float32
}{
	{SeedSet: true,
		SetSeed:    5,
		RandSource: rand.New(rand.NewSource(10)),
		Expected:   Float32Base{A: 0.24291174, C: 0.17921568, G: 0.39697546, T: 0.18089716},
		Precision:  1e-3,
	}, {SeedSet: false,
		SetSeed:    3,
		RandSource: nil,
		Expected:   Float32Base{A: 0.23355502, C: 0.21170676, G: 0.30556238, T: 0.24917582},
		Precision:  1e-3,
	},
}

func TestRandBase(t *testing.T) {
	var observed Float32Base
	for _, v := range RandBaseTests {
		observed = RandBase(v.SeedSet, v.SetSeed, v.RandSource)
		if !EqualBase(observed, v.Expected, v.Precision) {
			t.Errorf("Error: pDna randBase not as expected.\n")
		}
	}
}
