package pDna

import (
	"testing"
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

var DiffTests = []struct {
	P        Float32Base
	Q        Float32Base
	Expected Float64Diff
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
		Expected: Float64Diff{
			A: 0,
			C: 0,
			G: 0,
			T: 0,
		},
	},
	{P: Float32Base{
		A: 0.3,
		C: 0.2,
		G: 0.5,
		T: 0,
	},
		Q: Float32Base{
			A: 0.5,
			C: 0.1,
			G: 0.3,
			T: 0.1,
		},
		Expected: Float64Diff{
			A: -0.19999998807907104,
			C: 0.10000000149011612,
			G: 0.19999998807907104,
			T: -0.10000000149011612,
		},
	},
}

func TestDiff(t *testing.T) {
	var observed Float64Diff
	for _, v := range DiffTests {
		observed = Diff(v.P, v.Q)
		if observed != v.Expected {
			t.Errorf("Error: in pDna, Diff. Output was not as expected. observed: %v, expected: %v\n", observed, v.Expected)
		}
	}
}

var MagTests = []struct {
	D        Float64Diff
	Expected float64
}{
	{D: Float64Diff{
		A: 0,
		C: 0,
		G: 0,
		T: 0,
	},
		Expected: 0,
	},
	{D: Float64Diff{
		A: -0.19999998807907104,
		C: 0.10000000149011612,
		G: 0.19999998807907104,
		T: -0.10000000149011612,
	},
		Expected: 0.31622775188035535,
	},
}

func TestMag(t *testing.T) {
	var observed float64
	for _, v := range MagTests {
		observed = Mag(v.D)
		if observed != v.Expected {
			t.Errorf("Error: in pDna, Mag. Output was not as expected. observed: %v, expected: %v\n", observed, v.Expected)
		}
	}
}

var DistTests = []struct {
	P        Float32Base
	Q        Float32Base
	Expected float64
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
		Expected: 0,
	},
	{P: Float32Base{
		A: 0.3,
		C: 0.2,
		G: 0.5,
		T: 0,
	},
		Q: Float32Base{
			A: 0.5,
			C: 0.1,
			G: 0.3,
			T: 0.1,
		},
		Expected: 0.31622775188035535,
	},
}

func TestDist(t *testing.T) {
	var observed float64
	for _, v := range DistTests {
		observed = Dist(v.P, v.Q)
		if observed != v.Expected {
			t.Errorf("Error: in pDna, Dist. Output was not as expected. observed: %v, expected: %v\n", observed, v.Expected)
		}
	}
}
