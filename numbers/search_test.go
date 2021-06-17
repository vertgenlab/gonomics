package numbers

import (
	"testing"
)

func negativeNormalClosure(mu float64, sigma float64) func(float64) float64 {
	return func(x float64) float64 {
		return -1 * NormalDist(x, mu, sigma)
	}
}

func negativeGammaClosure(alpha float64, beta float64) func(float64) float64 {
	return func(x float64) float64 {
		return -1 * GammaDist(x, alpha, beta)
	}
}

var GgsMaxTests = []struct {
	f func(float64) float64
	a float64
	b float64
	epsilon float64
	expected float64
}{
	{NormalClosure(1.0, 0.1), -10, 10, 1e-5, 1.0},
	{GammaClosure(7.0, 3.0), 0, 20, 1e-5, 2.0},
}

var GgsMinTests = []struct {
	f func(float64) float64
	a float64
	b float64
	epsilon float64
	expected float64
}{
	{negativeNormalClosure(1.0, 0.1), -10.0, 10.0, 1e-5, 1.0},
	{negativeGammaClosure(7.0, 3.0), 0, 20, 1e-5, 2.0},
}

func TestGoldenSectionMinSearch(t *testing.T) {
	for _, test := range GgsMinTests {
		actual := GoldenSectionMinSearch(test.f, test.a, test.b, test.epsilon)
		if actual > test.expected*1.01 || actual < test.expected*0.99 {
			t.Errorf("Error in TestGoldenSectionMinSearch: expected %e. actual: %e.", test.expected, actual)
		}
	}
}

func TestGoldenSectionMaxSearch(t *testing.T) {
	for _, test := range GgsMaxTests {
		actual := GoldenSectionMaxSearch(test.f, test.a, test.b, test.epsilon)
		if actual > test.expected*1.01 || actual < test.expected*0.99 {
			t.Errorf("Error in TestGoldenSectionMaxSearch: expected %e. actual: %e.", test.expected, actual)
		}
	}
}
