package fit

import "testing"

var LagrangeInterpolationTests = []struct {
	queryX     float64
	dataPoints [][]float64
	expectedY  float64
}{
	{ //first parabola is y=x^2 + 1
		queryX: 0.5,
		dataPoints: [][]float64{
			{1.0, 2.0},
			{2.0, 5.0},
			{3.0, 10.0},
		},
		expectedY: 1.25,
	},
	{ // next is linear interpolation y= 3x + 2
		queryX: -1.0,
		dataPoints: [][]float64{
			{4.0, 14.0},
			{2.0, 8.0},
		},
		expectedY: -1.0,
	},
	{ // testing a cubic function y = 2x^3 + 6x^2 + 4x + 3
		queryX: 5.0,
		dataPoints: [][]float64{
			{-2.4311, 0.0}, //the root
			{-1.577, 3.77}, // local maximum
			{-0.423, 2.23}, // local minimum
			{0, 3},         // intercept
		},
		expectedY: 423,
	},
}

func TestLagrangeInterpolation(t *testing.T) {
	var currY float64
	for _, v := range LagrangeInterpolationTests {
		currY = LagrangeInterpolation(v.queryX, v.dataPoints)
		if (currY-v.expectedY)/v.expectedY > 0.05 {
			t.Errorf("Error: interpolated value: %v is not as expected: %v.\n", currY, v.expectedY)
		}
	}
}
