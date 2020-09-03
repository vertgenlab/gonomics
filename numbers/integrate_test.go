package numbers

import (
	"math"
	"testing"
	//DEBUG: "fmt"
)

var definiteIntegralTests = []struct {
	f      func(float64) float64
	a      float64
	b      float64
	answer float64
}{
	{func(x float64) float64 { return x }, 0, 1, 0.5},
	{func(x float64) float64 { return math.Pow(x, 2) + 1 }, 0, 2, 14.0 / 3.0},
	{func(x float64) float64 { return 6*math.Pow(x, 2) - 5*x + 2 }, -3, 1, 84},
	{func(x float64) float64 { return 1 / (math.Sqrt(math.Pi)) * math.Exp(-1*math.Pow(x, 2)) }, 0, 1, 0.42135039647485743467},
}

func TestDefiniteIntegral(t *testing.T) {
	for _, test := range definiteIntegralTests {
		calculated := DefiniteIntegral(test.f, test.a, test.b)
		if calculated > test.answer*1.01 || calculated < test.answer*0.99 {
			t.Errorf("For a definite integral we expected %e, but got %e", test.answer, calculated)
		}
	}
}

func BenchmarkDefiniteIntegral(b *testing.B) {
	for n := 0; n < b.N; n++ {
		for _, test := range definiteIntegralTests {
			DefiniteIntegral(test.f, test.a, test.b)
		}
	}
}

//answer from WolframAlpha
var logIntegralTests = []struct {
	f      func(float64) float64
	a      float64
	b      float64
	answer float64
}{
//	{func(x float64) float64 { return x }, 2, 3, 2.62},
	{func(x float64) float64 { return x }, 2, 5, 4.94893},
	{func(x float64) float64 { return math.Log(x*x) }, 2, 11, 6.08904},
}

func TestLogIntegral(t *testing.T) {
	for _, test := range logIntegralTests {
		calculated := LogIntegrate(test.f, test.a, test.b, 1000)
		//DEBUG: fmt.Printf("Calculated: %f. \n", calculated)
		if calculated > test.answer*1.01 || calculated < test.answer*0.99 {
			t.Errorf("For a logSpace integral we expected %e, but got %e", test.answer, calculated)
		}
	}
}

func TestLogIntegrateIterative(t *testing.T) {
	for _, test := range logIntegralTests {
		calculated := LogIntegrateIterative(test.f, test.a, test.b, 30, 1e-8)
		if calculated > test.answer*1.00001 || calculated < test.answer*0.99999 {
				t.Errorf("For a logSpace integral we expected %e, but got %e", test.answer, calculated)
			}
	}
}
