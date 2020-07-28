package numbers

import (
	"testing"
)

//pulls 100000 values from a gamma distribution and passes if the mean and variance of the values is within 10 percent of the expected values.
func TestRandGamma(t *testing.T) {
	alpha1 := 1.0
	beta1 := 1.0
	alpha2 := 5.0
	beta2 := 1.0
	alpha3 := 2.0
	beta3 := 6.0

	expectedMean1 := 1.0
	expectedMean2 := 5.0
	expectedMean3 := 1.0 / 3.0
	expectedVariance1 := 1.0
	expectedVariance2 := 5.0
	expectedVariance3 := 1.0 / 18.0

	var list1, list2, list3 []float64
	list1 = make([]float64, 100000)
	list2 = make([]float64, 100000)
	list3 = make([]float64, 100000)

	for i := 0; i < 100000; i++ {
		list1[i] = RandGamma(alpha1, beta1)
		list2[i] = RandGamma(alpha2, beta2)
		list3[i] = RandGamma(alpha3, beta3)
	}
	ave1 := AverageFloat64(list1)
	var1 := VarianceFloat64(list1)
	ave2 := AverageFloat64(list2)
	var2 := VarianceFloat64(list2)
	ave3 := AverageFloat64(list3)
	var3 := VarianceFloat64(list3)

	if ave1 < expectedMean1 * 0.9 || ave1 > expectedMean1 * 1.1 {
		t.Errorf("Gamma distribution simulated average outside acceptable range. Input : %e. Expected: %e.", ave1, expectedMean1)
	}
	if ave2 < expectedMean2 * 0.9 || ave2 > expectedMean2 * 1.1 {
		t.Errorf("Gamma distribution simulated average outside acceptable range. Input : %e. Expected: %e.", ave2, expectedMean2)
	}
	if ave3 < expectedMean3 * 0.9 || ave3 > expectedMean3 * 1.1 {
		t.Errorf("Gamma distribution simulated average outside acceptable range. Input : %f. Expected: %f. Accepted range: %f -- %f", ave3, expectedMean3, expectedMean3 * 0.9, expectedMean3 * 1.1)
	}
	if var1 < expectedVariance1 * 0.9 || var1 > expectedVariance1 * 1.1 {
		t.Errorf("Gamma distribution simulated variance outside acceptable range. Input : %e. Expected: %e.", var1, expectedVariance1)
	}
	if var2 < expectedVariance2 * 0.9 || var2 > expectedVariance2 * 1.1 {
		t.Errorf("Gamma distribution simulated variance outside acceptable range. Input : %e. Expected: %e.", var2, expectedVariance2)
	}
	if var3 < expectedVariance3 * 0.9 || var3 > expectedVariance3 * 1.1 {
		t.Errorf("Gamma distribution simulated variance outside acceptable range. Input : %e. Expected: %e.", var3, expectedVariance3)
	}
}