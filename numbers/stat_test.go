package numbers

import (
	"math"
	"testing"
)

var Float64Tests = []struct {
	Input    []float64
	Average  float64
	Variance float64
	StdDev   float64
}{
	{
		Input:    []float64{5, 8, 8, 3, 9, 7, 6, 7, 12, 9, 3, 1},
		Average:  6.5,
		Variance: 9.545454545454545,
		StdDev:   3.089571903266623,
	},
}

func TestFloat64(t *testing.T) {
	var ave, variance, stdDev float64
	for _, v := range Float64Tests {
		ave = AverageFloat64(v.Input)
		if ave != v.Average {
			t.Errorf("Error in AverageFloat64. Expected: %v. Output: %v.", v.Average, ave)
		}
		variance = VarianceFloat64(v.Input)
		if variance != v.Variance {
			t.Errorf("Error in VarianceFloat64. Expected: %v. Output: %v.", v.Variance, variance)
		}
		stdDev = StandardDeviationFloat64(v.Input)
		if stdDev != v.StdDev {
			t.Errorf("Error in StandardDeviationFloat64. Expected: %v. Output: %v.", v.StdDev, stdDev)
		}
	}
}

var PearsonTests = []struct {
	a        []float64
	b        []float64
	expected float64
}{
	{
		a:        []float64{3, 4, 5, 6, 7, 7, 6, 5, 4, 3, 5},
		b:        []float64{3, 4, 5, 6, 7, 7, 6, 5, 4, 3, 5},
		expected: 1,
	},
	{
		a:        []float64{3, 2, 5, 12, 7, 7, 6, 3, 4, 3, 5},
		b:        []float64{3, 4, 5, 6, 7, 7, 6, 5, 4, 3, 5},
		expected: 0.7015963532704358,
	},
}

func TestPearson(t *testing.T) {
	var observed float64
	for _, v := range PearsonTests {
		observed = Pearson(v.a, v.b)
		if math.Abs(observed-v.expected)/v.expected > 0.0001 {
			t.Errorf("Error in Pearson. Expected: %v. Observed: %v.", v.expected, observed)
		}
	}
}

func TestBenjaminiHochberg(t *testing.T) {
	var fail bool = false
	var input [][]float64 = [][]float64{{0, 0.005}, {1, 0.009}, {2, 0.019}, {3, 0.022}, {4, 0.051}, {5, 0.101}, {6, 0.361}, {7, 0.387}}
	var exp [][]float64 = [][]float64{{0, 0.036}, {1, 0.036}, {2, 0.044}, {3, 0.044}, {4, 0.0816}, {5, 0.135}, {6, 0.387}, {7, 0.387}}
	out := BenjaminiHochberg(input)
	for i := range out {
		if exp[i][0] != out[i][0] {
			fail = true
		}
		if math.Abs(exp[i][1]-out[i][1]) > 0.01 {
			fail = true
		}
	}
	if fail {
		t.Errorf("error in BenjaminiHochberg. Expected pValues: %v\nObserved pValues: %v", exp, out)
	}
}
