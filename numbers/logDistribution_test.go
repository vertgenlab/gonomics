package numbers

import (
	"math"
	"testing"
)

var BinomialDistLogTests = []struct {
	n      int
	k      int
	p      float64
	answer float64
}{
	{20, 4, 0.6, -8.21825},
	{20, 20, 0.6, -10.2165},
}

func TestBinomialDistLog(t *testing.T) {
	for _, test := range BinomialDistLogTests {
		calculated := BinomialDistLog(test.n, test.k, test.p)
		if math.Abs(calculated-test.answer)/test.answer > 0.000001 {
			t.Errorf("BinomialDistLog for n: %v, k: %v, and p: %f returned %f. Expected: %f.", test.n, test.k, test.p, calculated, test.answer)
		}
	}
}

var AddBinomMapEntryTests = []struct {
	n int
}{
	{4},
	{6},
}

func TestAddBinomMapEntry(t *testing.T) {
	for _, test := range AddBinomMapEntryTests {
		calculated := AddBinomMapEntry(test.n)
		for k := 0; k < len(calculated); k++ {
			if calculated[k] != BinomCoefficientLog(test.n, k) {
				t.Errorf("AddBinomMapEntry for n: %v, k: %v returned %f. Expected: %f.", test.n, k, calculated[k], BinomCoefficientLog(test.n, k))
			}
		}
	}
}

var BinomialDistLogMapTests = []struct {
	n      int
	k      int
	p      float64
	answer float64
}{
	{20, 4, 0.6, -8.21825},
	{20, 20, 0.6, -10.216512},
	{30, 14, 0.6, -3.01706},
}

//tests if all of the values are correct and tests the number of entries in the binomMap at the end
func TestBinomalDistLogMap(t *testing.T) {
	binomMap := make(map[int][]float64)
	for _, test := range BinomialDistLogMapTests {
		calculated := BinomialDistLogMap(test.n, test.k, test.p, &binomMap)
		if (calculated-test.answer)/test.answer > 0.000001 {
			t.Errorf("BinomialDistLogMap for n: %v, k: %v, and p: %f returned %f. Expected: %f.", test.n, test.k, test.p, calculated, test.answer)
		}
	}
	if len(binomMap) != 2 {
		t.Errorf("BinomialDistLogMap had binomMap of length %v. Expected 2.", len(binomMap))
	}
}

/*
var binomMap map[int][]float64
var BinomialDistLogMapTests = []struct {

}*/
