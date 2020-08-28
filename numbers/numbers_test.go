package numbers

import (
	"fmt"
	"math"
	"testing"
)

var fisherExactLessTests = []struct {
	a      int     // top left in matrix
	b      int     // top right in matrix
	c      int     // bottom left in matrix
	d      int     // bottom right in matrix
	pvalue float64 // pvalue of a being less
}{
	{0, 3, 7, 10, 0.2509},
	{1, 3, 9, 10, 0.4037},
	{1, 70, 13, 800, 0.6882},
	{2, 70, 13, 800, 0.8843},
	{3, 70, 13, 800, 0.9639},
	{3, 300, 700, 8000, 4.487e-8},
	{90, 101, 700, 8000, 1},
}

var fisherExactGreaterTests = []struct {
	a      int     // top left in matrix
	b      int     // top right in matrix
	c      int     // bottom left in matrix
	d      int     // bottom right in matrix
	pvalue float64 // pvalue of a being greater
}{
	{73, 2076, 5057, 180930, 0.0353347421506263},
	{122, 3709, 5008, 179297, 0.0463757733199719},
	{89, 2777, 5041, 180229, 0.117081148851346},
	{105, 3137, 5025, 179869, 0.0427107523491912},
	{86, 2962, 5044, 180044, 0.388336598185756},
	{18, 297, 5112, 182709, 0.00287720449427838},
	{20, 328, 5110, 182678, 0.00163411525648603},
	{107, 3046, 5023, 179960, 0.0138585128159474},
	{108, 3297, 5022, 179709, 0.0623694515380971},
	{59, 2476, 5071, 180530, 0.906449513485591},
	{111, 3357, 5019, 179649, 0.0494041701774109},
	{0, 0, 94, 90, 1},
	{108, 1432, 742, 70208, 2.95e-50},
	{76, 542, 774, 71098, 1.04e-52},
	{47, 1253, 50, 84636, 7.12e-59},
	{112, 2417, 87, 114989, 2.40e-131},
	{313, 1977, 537, 69663, 6.59e-245},
	{520, 1239, 889, 97088, math.SmallestNonzeroFloat64},
}

func TestFisherLess(t *testing.T) {
	for _, test := range fisherExactLessTests {
		calculated := FisherExact(test.a, test.b, test.c, test.d, true)
		if calculated > test.pvalue*1.01 || calculated < test.pvalue*0.99 {
			t.Errorf("For a fisher test (less) on (%d, %d, %d, %d): expected %e, but got %e", test.a, test.b, test.c, test.d, test.pvalue, calculated)
		}
	}
}

func TestFisherGreater(t *testing.T) {
	for _, test := range fisherExactGreaterTests {
		calculated := FisherExact(test.a, test.b, test.c, test.d, false)
		if calculated > test.pvalue*1.01 || calculated < test.pvalue*0.99 {
			t.Errorf("For a fisher test (greater) on (%d, %d, %d, %d): expected %e, but got %e", test.a, test.b, test.c, test.d, test.pvalue, calculated)
		}
	}
}

var BinomCoefficientTests = []struct {
	n      int
	k      int
	answer int
}{
	{0, 0, 1},
	{9, 0, 1},
	{3, 2, 3},
	{30, 12, 86493225},
	{13, 13, 1},
	{9, 2, 36},
	{11, 5, 462},
}

//answers from WolframAlpha, which can handle large choose inputs
var BinomCoefficientLogTests = []struct {
	n      int
	k      int
	answer float64
}{
	{2000, 1000, 1382.267994},
	{10000, 8000, 4999.416374},
}

func TestBinomCoefficient(t *testing.T) {
	for _, test := range BinomCoefficientTests {
		calculated := BinomCoefficient(test.n, test.k)
		if calculated != test.answer {
			t.Errorf("For BinomialCoefficient(%d, %d) we would expect %d, but we got %d", test.n, test.k, test.answer, calculated)
		}
	}
}

func TestBinomCoefficientLog(t *testing.T) {
	//first we see if we have consistency with the non-log version
	for _, test := range BinomCoefficientTests {
		calculated := math.Exp(BinomCoefficientLog(test.n, test.k))
		if fmt.Sprintf("%f", calculated) != fmt.Sprintf("%f", float64(test.answer)) {
			t.Errorf("ForBinomialCoefficientLog(%d, %d) we would expect %f, but we got %f", test.n, test.k, float64(test.answer), calculated)
		}
	}
	//now we test big numbers
	for _, test := range BinomCoefficientLogTests {
		calculated := BinomCoefficientLog(test.n, test.k)
		if fmt.Sprintf("%f", calculated) != fmt.Sprintf("%f", test.answer) {
			t.Errorf("ForBinomialCoefficientLog(%d, %d) we would expected %f, but we got %f", test.n, test.k, test.answer, calculated)
		}
	}
}
