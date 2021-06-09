package numbers

import (
	"math"
	"testing"
)

var highestDensityIntervalTests = []struct {
	Trace McmcTrace
	Start float64
	End float64
	Mean float64
	Proportion float64
}{
		{McmcTrace{Parameter: []float64{0.1, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 20.0}}, 2.0, 9.0, 6.41,0.8},
}

func TestHighestDensityInterval(t *testing.T) {
	for _, test := range highestDensityIntervalTests {
		s, e := HighestDensityInterval(test.Trace, test.Proportion)
		if s != test.Start || e != test.End {
			t.Errorf("Error in HighestDensityInterval. Start: %f. End: %f. ExpectedStart:%f. ExpectedEnd:%f.", s, e, test.Start, test.End)
		}
	}
}

func TestMeanMcmcTrace(t *testing.T) {
	for _, test := range highestDensityIntervalTests {
		m := MeanMcmcTrace(test.Trace)
		if math.Abs(m - test.Mean) > 0.1 {
			t.Errorf("Error in MeanMcmcTrace. Mean: %f. Expected Mean: %f.", m, test.Mean)
		}
	}
}