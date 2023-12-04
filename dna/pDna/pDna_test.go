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
