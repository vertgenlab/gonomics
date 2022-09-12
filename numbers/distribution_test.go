package numbers

import (
	"fmt"
	"math"
	"testing"
)

//All expected values are calculated from online calculator tools or in the R programming language.

func ExplicitBinomialDist(n int, k int, p float64) float64 {
	return float64(BinomCoefficient(n, k)) * math.Pow(p, float64(k)) * math.Pow(1-p, float64(n-k))
}

func TestBinomialDist(t *testing.T) {
	input1 := ExplicitBinomialDist(20, 4, 0.6)
	expected1, _ := BinomialDist(20, 4, 0.6)
	input2 := ExplicitBinomialDist(20, 20, 0.6)
	expected2, _ := BinomialDist(20, 20, 0.6)
	input3 := ExplicitBinomialDist(20, 0, 0.6)
	expected3, _ := BinomialDist(20, 0, 0.6)
	if fmt.Sprintf("%e", input1) != fmt.Sprintf("%e", expected1) {
		t.Errorf("Do not match. Input : %e. Expected: %e.", input1, expected1)
	}
	if fmt.Sprintf("%e", input2) != fmt.Sprintf("%e", expected2) {
		t.Errorf("Do not match. Input : %e. Expected: %e.", input2, expected2)
	}
	if fmt.Sprintf("%e", input3) != fmt.Sprintf("%e", expected3) {
		t.Errorf("Do not match. Input : %e. Expected: %e.", input3, expected3)
	}
}

func TestBinomialSum(t *testing.T) {
	input1 := BinomialLeftSummation(20, 1, 0.6)
	expected1 := 3.408486e-07
	input2 := BinomialLeftSummation(20, 20, 0.6)
	expected2 := 1.000000e+00
	input3 := BinomialRightSummation(20, 4, 0.6)
	expected3 := 9.999527e-01
	input4 := BinomialRightSummation(20, 16, 0.4)
	expected4 := 3.170311e-04
	if fmt.Sprintf("%e", input1) != fmt.Sprintf("%e", expected1) {
		t.Errorf("Do not match. Input : %e. Expected: %e.", input1, expected1)
	}
	if fmt.Sprintf("%e", input2) != fmt.Sprintf("%e", expected2) {
		t.Errorf("Do not match. Input : %e. Expected: %e.", input2, expected2)
	}
	if fmt.Sprintf("%e", input3) != fmt.Sprintf("%e", expected3) {
		t.Errorf("Do not match. Input : %e. Expected: %e.", input3, expected3)
	}
	if fmt.Sprintf("%e", input4) != fmt.Sprintf("%e", expected4) {
		t.Errorf("Do not match. Input : %e. Expected: %e.", input4, expected4)
	}
}

func TestPoissonDist(t *testing.T) {
	input1 := PoissonDist(4, 5)
	expected1 := 1.754674e-01
	input2 := PoissonDist(0, 5)
	expected2 := 6.737947e-03
	if fmt.Sprintf("%e", input1) != fmt.Sprintf("%e", expected1) {
		t.Errorf("Do not match. Input : %e. Expected: %e.", input1, expected1)
	}
	if fmt.Sprintf("%e", input2) != fmt.Sprintf("%e", expected2) {
		t.Errorf("Do not match. Input : %e. Expected: %e.", input2, expected2)
	}
}

func TestPoissonSum(t *testing.T) {
	input1 := PoissonLeftSummation(4, 5)
	expected1 := 4.404933e-01
	input2 := PoissonLeftSummation(0, 5)
	expected2 := 6.737947e-03
	input3 := PoissonRightSummation(7, 5)
	expected3 := 2.378165e-01
	input4 := PoissonRightSummation(0, 5)
	expected4 := 1.000000e+00
	if fmt.Sprintf("%e", input1) != fmt.Sprintf("%e", expected1) {
		t.Errorf("Do not match. Input : %e. Expected: %e.", input1, expected1)
	}
	if fmt.Sprintf("%e", input2) != fmt.Sprintf("%e", expected2) {
		t.Errorf("Do not match. Input : %e. Expected: %e.", input2, expected2)
	}
	if fmt.Sprintf("%e", input3) != fmt.Sprintf("%e", expected3) {
		t.Errorf("Do not match. Input : %e. Expected: %e.", input3, expected3)
	}
	if fmt.Sprintf("%e", input4) != fmt.Sprintf("%e", expected4) {
		t.Errorf("Do not match. Input : %e. Expected: %e.", input4, expected4)
	}
}

func TestNormalDist(t *testing.T) {
	input := NormalDist(0, 0, 1)
	expected := 3.989423e-01
	if fmt.Sprintf("%e", input) != fmt.Sprintf("%e", expected) {
		t.Errorf("Do not match. Input : %e. Expected: %e.", input, expected)
	}
}

//TODO: Romberg's method fails when evaluating large probability regions.
func TestNormalIntegral(t *testing.T) {
	f := NormalClosure(0, 1)
	input1 := DefiniteIntegral(f, 3, 200)
	expected1 := 1.349890e-03
	input2 := DefiniteIntegral(f, -200, -6)
	expected2 := 6.389367e-10
	if fmt.Sprintf("%e", input1) != fmt.Sprintf("%e", expected1) {
		t.Errorf("Do not match. Input : %e. Expected: %e.", input1, expected1)
	}
	if fmt.Sprintf("%e", input2) != fmt.Sprintf("%e", expected2) {
		t.Errorf("Do not match. Input : %e. Expected: %e.", input2, expected2)
	}
}

func TestBetaDist(t *testing.T) {
	input1 := BetaDist(0.3, 2, 3)
	expected1 := 1.764000e+00
	input2 := BetaDist(0, 2, 3)
	expected2 := 0.000000e+00
	if fmt.Sprintf("%e", input1) != fmt.Sprintf("%e", expected1) {
		t.Errorf("Do not match. Input : %e. Expected: %e.", input1, expected1)
	}
	if fmt.Sprintf("%e", input2) != fmt.Sprintf("%e", expected2) {
		t.Errorf("Do not match. Input : %e. Expected: %e.", input2, expected2)
	}
}

func TestBetaIntegral(t *testing.T) {
	input1 := BetaRightIntegral(0.3, 2, 3)
	expected1 := 6.517000e-01
	input2 := BetaLeftIntegral(0.3, 2, 3)
	expected2 := 3.483000e-01
	if fmt.Sprintf("%e", input1) != fmt.Sprintf("%e", expected1) {
		t.Errorf("Do not match. Input : %e. Expected: %e.", input1, expected1)
	}
	if fmt.Sprintf("%e", input2) != fmt.Sprintf("%e", expected2) {
		t.Errorf("Do not match. Input : %e. Expected: %e.", input2, expected2)
	}
}

func TestGammaDist(t *testing.T) {
	input1 := GammaDist(0.2, 4, 10)
	expected1 := 1.804470e+00
	input2 := GammaDist(3.5, 4, 2)
	expected2 := 1.042585e-01
	if fmt.Sprintf("%e", input1) != fmt.Sprintf("%e", expected1) {
		t.Errorf("Do not match. Input : %e. Expected: %e.", input1, expected1)
	}
	if fmt.Sprintf("%e", input2) != fmt.Sprintf("%e", expected2) {
		t.Errorf("Do not match. Input : %e. Expected: %e.", input2, expected2)
	}
}

//TODO: Similar romberg's problems: integrating over large ranges leads to failure of convergence.
func TestGammaIntegral(t *testing.T) {
	input1 := GammaLeftIntegral(3, 4, 2)
	expected1 := 8.487961e-01
	input2 := GammaRightIntegral(8, 4, 2)
	expected2 := 9.315161e-05
	if fmt.Sprintf("%e", input1) != fmt.Sprintf("%e", expected1) {
		t.Errorf("Do not match. Input : %e. Expected: %e.", input1, expected1)
	}
	if fmt.Sprintf("%e", input2) != fmt.Sprintf("%e", expected2) {
		t.Errorf("Do not match. Input : %e. Expected: %e.", input2, expected2)
	}
}

//TODO: Testing for KL Divergence is incomplete
//Expected values are calculated with the R package LaplacesDemon using the KLD function.
func TestKullbackLeiblerDivergence(t *testing.T) {
	//input1 := ContinuousKullbackLeiblerDivergence(NormalClosure(0, 1), NormalClosure(1, 1), -100, 100)
	//fmt.Printf("%v\n", input1)
}
