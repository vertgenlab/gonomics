package fit

import "log"

// LagrangeInterpolation returns f(queryX) for a given queryX. The function f(x) is
// estimated by Lagrange Interpolation using a set of input dataPoints in [][]float64 format,
// each data point is a []float64 with two elements, the first value is the x value and the
// second value is the y value.
// The estimated function is of order len(dataPoints)-1. in other words, two input points will
// define a linear interpolation, three input points will define a quadratic interpolation, and so on.
// Be wary of extrapolation.
func LagrangeInterpolation(queryX float64, dataPoints [][]float64) float64 {
	var answer float64 = 0
	var currProduct float64

	if len(dataPoints) < 2 {
		log.Fatalf("Error: require at least two data points to perform Lagrange Interpolation.")
	}

	for currDataPoint := range dataPoints {
		currProduct = 1.0
		for currOtherPoint := range dataPoints {
			if currDataPoint != currOtherPoint {
				currProduct *= (queryX - dataPoints[currOtherPoint][0]) / (dataPoints[currDataPoint][0] - dataPoints[currOtherPoint][0])
			}
		}
		answer += dataPoints[currDataPoint][1] * currProduct
	}
	return answer
}
