package numbers

import (
	"fmt"
	"math"
	"math/rand"
)

const logProbEpsilon = 1e-08

type MixtureModel struct {
	Data              []float64 // 1d data slice
	K                 int       // number of component distributions
	Means             []float64 // means for each component. len(means) == k
	Variances         []float64 // variances for each component. len(variances) == k
	Weights           []float64
	MaxIter           int         // maximum number of iterations for EM step. 0 is until convergence
	NLogLikelihood    float64     // negative likelihood to be minimized
	responsibility    [][]float64 // first index is component, second index is data point
	responsibilitySum []float64   // sum of responsibilities above. len(responsibilitySum) == k
}

func initMixtureModel(data []float64, k int, maxIterations int, mm *MixtureModel) {
	if mm == nil {
		mm = new(MixtureModel)
	}
	mm.MaxIter = maxIterations
	mm.Data = data
	mm.K = k
	mm.NLogLikelihood = math.MaxFloat64
	if cap(mm.Means) >= k {
		mm.Means = mm.Means[0:k]
	} else {
		mm.Means = make([]float64, k)
	}
	if cap(mm.Variances) >= k {
		mm.Variances = mm.Variances[0:k]
	} else {
		mm.Variances = make([]float64, k)
	}

	// TODO smarter initial guess for mean and variance (k-means/PCA)
	for i := range mm.Means {
		mm.Means[i] = rand.Float64() * 100
		mm.Variances[i] = 100
	}

	if cap(mm.responsibility) >= k {
		mm.responsibility = mm.responsibility[0:k]
	} else {
		mm.responsibility = make([][]float64, k)
	}

	for i := range mm.responsibility {
		if cap(mm.responsibility[i]) >= len(data) {
			mm.responsibility[i] = mm.responsibility[i][0:len(data)]
		} else {
			mm.responsibility[i] = make([]float64, len(data))
		}
	}

	if cap(mm.responsibilitySum) >= k {
		mm.responsibilitySum = mm.responsibilitySum[0:k]
	} else {
		mm.responsibilitySum = make([]float64, k)
	}

	if cap(mm.Weights) >= k {
		mm.Weights = mm.Weights[0:k]
	} else {
		mm.Weights = make([]float64, k)
	}

	for i := range mm.Weights {
		mm.Weights[i] = 1 / float64(k)
	}
}

func RunMixtureModel(data []float64, k int, maxIterations int, mm *MixtureModel) {
	initMixtureModel(data, k, maxIterations, mm)

	var converged bool
	for iter := 0; !converged && iter < mm.MaxIter; iter++ {

		// E step
		var model float64
		for i := range mm.Data {
			var totalProb float64
			for j := 0; j < mm.K; j++ {
				var prob float64
				prob = mm.Weights[j] * calculate1dProb(mm.Data[i], mm.Means[j], mm.Variances[j])
				totalProb += prob
				mm.responsibility[j][i] = prob
			}

			if totalProb > 0 { // normalize responsibilities for each datapoint
				for j := 0; j < mm.K; j++ {
					mm.responsibility[j][i] /= totalProb
				}
			}
			model -= math.Log(totalProb)
		}

		// Check that model negative log likelihood is decreasing
		// (negative log likelihood is guaranteed to be non-increasing).
		// If the current value does not differ from the previous,
		// the algorithm has converged (possibly to a local minimum).
		//fmt.Println("Iteration number:", iter)
		//fmt.Println(model)
		//fmt.Println(mm.NLogLikelihood)
		//fmt.Println(mm.Variances)
		//fmt.Println(mm.Means)
		if mm.NLogLikelihood-model < logProbEpsilon {
			converged = true
			fmt.Println("Converged on iteration:", iter)
		}
		mm.NLogLikelihood = model

		// M step
		// sum responsibilities of each data point to each component
		resetResSum(mm)
		for i := range mm.Data {
			for j := 0; j < mm.K; j++ {
				mm.responsibilitySum[j] += mm.responsibility[j][i]
			}
		}

		// normalize weights to 0-1
		for j := 0; j < mm.K; j++ {
			mm.Weights[j] = mm.responsibilitySum[j] / float64(len(mm.Data))
		}

		for j := 0; j < mm.K; j++ {
			mm.Means[j] = 0 // we will calculate a new one
			for i := range mm.Data {
				mm.Means[j] += mm.responsibility[j][i] * mm.Data[i]
			}

			if mm.responsibilitySum[j] > 0 {
				mm.Means[j] /= mm.responsibilitySum[j]
			}

			for i := range mm.Data {
				var diff float64
				diff = mm.Data[i] - mm.Means[j]
				mm.Variances[j] += mm.responsibility[j][i] * diff * diff
			}

			if mm.responsibilitySum[j] > 0 {
				mm.Variances[j] /= mm.responsibilitySum[j]
			}

			mm.Variances[j] = math.Sqrt(mm.Variances[j])
		}
	}
}

func calculate1dProb(val, mean, std float64) float64 {
	var frac, power float64
	frac = 1 / (std * math.Sqrt(2*math.Pi))
	power = -0.5 * math.Pow((val-mean)/std, 2)
	return frac * math.Exp(power)
}

func resetResSum(mm *MixtureModel) {
	for i := range mm.responsibilitySum {
		mm.responsibilitySum[i] = 0
	}
}
