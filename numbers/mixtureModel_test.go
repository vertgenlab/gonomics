package numbers

import (
	"fmt"
	"math"
	"math/rand"
	"testing"

	"golang.org/x/exp/slices"
)

func TestRunMixtureModel(t *testing.T) {
	maxIterations := 200
	maxResets := 10
	countPerComponent := 300
	acceptableError := 0.05
	var meanA, meanB, stdA, stdB float64
	mm := new(MixtureModel)
	var converged bool
	data := make([]float64, 0, countPerComponent*2)
	for i := 0; i < 100; i++ {
		meanA = (rand.Float64() * 100) + 20        // min mean of 20
		meanB = meanA + 20 + (rand.Float64() * 20) // at least 20 apart
		stdA = rand.Float64() * 10
		stdB = rand.Float64() * 10
		data = append(generateData(countPerComponent, meanA, stdA), generateData(countPerComponent, meanB, stdB)...)
		converged, _ = RunMixtureModel1D(data, 2, maxIterations, maxResets, 1e-10, mm)
		if !converged {
			continue
		}
		slices.Sort(mm.Means)
		if math.Abs(1-mm.Means[0]/meanA) > acceptableError || math.Abs(1-mm.Means[1]/meanB) > acceptableError {
			t.Errorf("problem with RunMixtureModel")
			fmt.Println(mm.Means, meanA, meanB, mm.Weights, converged)
		}
	}
}

func BenchmarkRunMixtureModel(b *testing.B) {
	countPerComponent := 100
	means := []float64{120, 480}
	std := []float64{20, 30}
	data := make([]float64, 0, 200)
	for i := range means {
		data = append(data, generateData(countPerComponent, means[i], std[i])...)
	}
	mm := new(MixtureModel)
	RunMixtureModel1D(data, 2, 200, 5, 1e-10, mm) // outside loop for initial memory allocation
	for i := 0; i < b.N; i++ {
		RunMixtureModel1D(data, 2, 200, 5, 1e-10, mm)
	}
}

// generate data from a normal distribution with noise.
func generateData(num int, mean, std float64) []float64 {
	ans := make([]float64, num)
	for i := range ans {
		ans[i] = (rand.NormFloat64() * std) + mean
	}
	return ans
}
