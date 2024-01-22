package wig

import (
	"fmt"
	"testing"
)

var PearsonTests = []struct {
	alpha             map[string]Wig
	beta              map[string]Wig
	missing           float64
	samplingFrequency float64
	expected          float64
}{
	{
		alpha: map[string]Wig{
			"chr1": {
				StepType: "fixedStep",
				Chrom:    "chr1",
				Start:    1,
				Step:     1,
				Values:   []float64{10, 5, 3, 5, 10},
			},
			"chr2": {
				StepType: "fixedStep",
				Chrom:    "chr2",
				Start:    1,
				Step:     1,
				Values:   []float64{2, -10, 4, 5, 6, 5, 4, 3, 7, 7, 9},
			},
		},
		beta: map[string]Wig{
			"chr2": {
				StepType: "fixedStep",
				Chrom:    "chr2",
				Start:    1,
				Step:     1,
				Values:   []float64{2, -10, 7, 5, 6, 5, 4, 3, 7, 7, 9},
			},
			"chr1": {
				StepType: "fixedStep",
				Chrom:    "chr1",
				Start:    1,
				Step:     1,
				Values:   []float64{10, 5, 3, 5, 10},
			},
		},
		missing:           -10,
		samplingFrequency: 1,
		expected:          0.951504453802689,
	},
}

func TestPearson(t *testing.T) {
	var observed float64
	for _, v := range PearsonTests {
		observed = Pearson(v.alpha, v.beta, v.missing, v.samplingFrequency)
		if fmt.Sprintf("%.6g", observed) != fmt.Sprintf("%.6g", v.expected) {
			t.Errorf("Error in wig.Pearson. Expected: %v. Output: %v.", v.expected, observed)
		}
	}
}
