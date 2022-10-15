package numbers

import (
	"math"
	"math/rand"
)

const logProbEpsilon = 1e-08

type MixtureModel struct {
	Data           []float64 // 1d data slice
	K              int       // number of component distributions
	Means          []float64 // means for each component. len(means) == k
	Stdev          []float64 // variances for each component. len(stdev) == k
	Weights        []float64
	MaxIter        int         // maximum number of iterations for EM step. 0 is until convergence
	LogLikelihood  float64     // negative likelihood to be minimized
	responsibility [][]float64 // first index is component, second index is data point
	posteriorsSum  []float64   // sum of responsibilities above. len(posteriorsSum) == k
	posteriors     [][]float64

	lamSigRatio    []float64
	logLamSigRatio []float64
	work           []float64
}

func initMixtureModel(data []float64, k int, maxIterations int, mm *MixtureModel) {
	if mm == nil {
		mm = new(MixtureModel)
	}
	mm.MaxIter = maxIterations
	mm.Data = data
	mm.K = k
	mm.LogLikelihood = math.MaxFloat64
	if cap(mm.Means) >= k {
		mm.Means = mm.Means[0:k]
	} else {
		mm.Means = make([]float64, k)
	}
	if cap(mm.Stdev) >= k {
		mm.Stdev = mm.Stdev[0:k]
	} else {
		mm.Stdev = make([]float64, k)
	}
	if cap(mm.lamSigRatio) >= k {
		mm.lamSigRatio = mm.lamSigRatio[0:k]
	} else {
		mm.lamSigRatio = make([]float64, k)
	}
	if cap(mm.logLamSigRatio) >= k {
		mm.logLamSigRatio = mm.logLamSigRatio[0:k]
	} else {
		mm.logLamSigRatio = make([]float64, k)
	}
	if cap(mm.work) >= k {
		mm.work = mm.work[0:k]
	} else {
		mm.work = make([]float64, k)
	}

	// TODO smarter initial guess for mean and variance (k-means/PCA)
	for i := range mm.Means {
		mm.Means[i] = rand.Float64() * 100
		mm.Stdev[i] = 10
	}
	mm.Means[0] = 0
	mm.Means[1] = 100

	if cap(mm.responsibility) >= k {
		mm.responsibility = mm.responsibility[0:k]
		mm.posteriors = mm.posteriors[0:k]
	} else {
		mm.responsibility = make([][]float64, k)
		mm.posteriors = make([][]float64, k)
	}

	for i := range mm.responsibility {
		if cap(mm.responsibility[i]) >= len(data) {
			mm.responsibility[i] = mm.responsibility[i][0:len(data)]
			mm.posteriors[i] = mm.posteriors[i][0:len(data)]
		} else {
			mm.responsibility[i] = make([]float64, len(data))
			mm.posteriors[i] = make([]float64, len(data))
		}
	}

	if cap(mm.posteriorsSum) >= k {
		mm.posteriorsSum = mm.posteriorsSum[0:k]
	} else {
		mm.posteriorsSum = make([]float64, k)
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
	var prevLogLikelihood float64
	for iter := 0; !converged && iter < mm.MaxIter; iter++ {

		// E step
		prevLogLikelihood = mm.LogLikelihood
		expectation(mm)
		//fmt.Println(mm.Means, mm.LogLikelihood, mm.Stdev)
		if iter > 2 && mm.LogLikelihood-prevLogLikelihood < logProbEpsilon {
			converged = true
			//fmt.Println("Converged on iteration:", iter)
		}

		// M step
		// sum responsibilities of each data point to each component
		maximization(mm)
	}
}

func resetResSum(mm *MixtureModel) {
	for i := range mm.posteriorsSum {
		mm.posteriorsSum[i] = 0
	}
}

// adapted from https://github.com/cran/mixtools/blob/master/src/normpost.c
func expectation(mm *MixtureModel) {
	var r, x, min, rowsum float64
	var i, j, minj int
	mm.LogLikelihood = -float64(len(mm.Data)/2) * 0.91893853320467274178 // -n/2 * log(2pi)
	for i = 0; i < mm.K; i++ {
		mm.lamSigRatio[i] = mm.Weights[i] / mm.Stdev[i]
		mm.logLamSigRatio[i] = math.Log(mm.lamSigRatio[i])
	}

	for i = range mm.Data {
		x = mm.Data[i]
		for j = 0; j < mm.K; j++ {
			r = x - mm.Means[j]
			r *= r
			mm.responsibility[j][i] = r
			r = r / (2 * mm.Stdev[j] * mm.Stdev[j])
			mm.work[j] = r

			/* Keep track of the smallest standardized squared residual.
			   By dividing everything by the component density with the
			   smallest such residual, the denominator of the posterior
			   is guaranteed to be at least one and cannot be infinite unless
			   the values of lambda or sigma are very large or small. This helps
			   prevent numerical problems when calculating the posteriors.*/
			if j == 0 || r < min {
				minj = j
				min = r
			}
		}
		/* At this stage, work contains the squared st'dized resids over 2 */
		rowsum = 1
		for j = 0; j < mm.K; j++ {
			if j == minj {
				mm.work[j] = 1
			} else {
				mm.work[j] = (mm.lamSigRatio[j] / mm.lamSigRatio[minj]) * math.Exp(min-mm.work[j])
				rowsum += mm.work[j]
			}
		}
		/* At this stage, work contains the normal density at data[i]
		   divided by the normal density with the largest st'dized resid
		   Thus, dividing by rowsum gives the posteriors: */
		for j = 0; j < mm.K; j++ {
			mm.posteriors[j][i] = mm.work[j] / rowsum
		}
		/* Finally, adjust the loglikelihood correctly */
		mm.LogLikelihood += math.Log(rowsum) - min + mm.logLamSigRatio[minj]
	}
}

func maximization(mm *MixtureModel) {
	resetResSum(mm)
	for i := range mm.Data {
		for j := 0; j < mm.K; j++ {
			mm.posteriorsSum[j] += mm.posteriors[j][i]
		}
	}

	// normalize weights to 0-1
	for j := 0; j < mm.K; j++ {
		mm.Weights[j] = mm.posteriorsSum[j] / float64(len(mm.Data))
	}

	for j := 0; j < mm.K; j++ {
		mm.Means[j] = 0
		for i := range mm.Data {
			mm.Means[j] += mm.posteriors[j][i] * mm.Data[i]
		}

		if mm.posteriorsSum[j] > 0 {
			mm.Means[j] /= mm.posteriorsSum[j]
		}

		for i := range mm.Data {
			mm.Stdev[j] += mm.posteriors[j][i] * mm.responsibility[j][i]
		}

		if mm.posteriorsSum[j] > 0 {
			mm.Stdev[j] /= mm.posteriorsSum[j]
		}

		mm.Stdev[j] = math.Sqrt(mm.Stdev[j])
	}
}

func oldExpectation(mm *MixtureModel) {
	var model float64
	for i := range mm.Data {
		var totalProb float64
		for j := 0; j < mm.K; j++ {
			var prob float64
			prob = mm.Weights[j] * calculate1dProb(mm.Data[i], mm.Means[j], mm.Stdev[j])
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
}

func oldMaximization(mm *MixtureModel) {
	resetResSum(mm)
	for i := range mm.Data {
		for j := 0; j < mm.K; j++ {
			mm.posteriorsSum[j] += mm.responsibility[j][i]
		}
	}

	// normalize weights to 0-1
	for j := 0; j < mm.K; j++ {
		mm.Weights[j] = mm.posteriorsSum[j] / float64(len(mm.Data))
	}

	for j := 0; j < mm.K; j++ {
		mm.Means[j] = 0 // we will calculate a new one
		for i := range mm.Data {
			mm.Means[j] += mm.responsibility[j][i] * mm.Data[i]
		}

		if mm.posteriorsSum[j] > 0 {
			mm.Means[j] /= mm.posteriorsSum[j]
		}

		for i := range mm.Data {
			var diff float64
			diff = mm.Data[i] - mm.Means[j]
			mm.Stdev[j] += mm.responsibility[j][i] * diff * diff
		}

		if mm.posteriorsSum[j] > 0 {
			mm.Stdev[j] /= mm.posteriorsSum[j]
		}

		mm.Stdev[j] = math.Sqrt(mm.Stdev[j])
	}
}

func calculate1dProb(val, mean, std float64) float64 {
	var frac, power float64
	frac = 1 / (std * math.Sqrt(2*math.Pi))
	power = -0.5 * math.Pow((val-mean)/std, 2)
	return frac * math.Exp(power)
}
