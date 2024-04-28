package pDna

import "testing"

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
