package numbers

import (
	"math"
)

// logProbEpsilon is the arbitrary set point for convergence. When the model likelihood is
// increasing by < logProbEpsilon, we say the model has converged, possibly to a local minimum.
const logProbEpsilon = 1e-08

// MixtureModel holds data, results, and working memory for running the EM algorithm.
type MixtureModel struct {
	Data           []float64   // 1d data slice
	K              int         // number of component distributions
	Means          []float64   // means for each component. len(means) == k
	Stdev          []float64   // variances for each component. len(stdev) == k
	Weights        []float64   // contribution of each gaussian to the model
	MaxIter        int         // maximum number of iterations for EM step. 0 is until convergence
	LogLikelihood  float64     // negative likelihood to be minimized
	residuals      [][]float64 // first index is component, second index is data point
	posteriors     [][]float64 // posterior values for each data point for each gaussian
	posteriorsSum  []float64   // sum of posteriors above. len(posteriorsSum) == k
	lamSigRatio    []float64   // holds the ratio of Weights to Stdev
	logLamSigRatio []float64   // holds the natural log of lamSigRatio
	work           []float64   // holds intermediate values for calculating LogLikelihood
}

// RunMixtureModel1D uses the expectation-maximization (EM) algorithm to find a mixture of k gaussian distributions that fit the input data slice.
// Note that this version of RunMixtureModel only works on 1d data. The EM algorithm works by iteratively refining the model until the performance
// of the model is no longer improving (i.e. it has converged). RunMixtureModel will iterate a maximum of maxIterations until retrying with new
// starting values until convergence or maxResets. RunMixtureModel will store the results of the model in mm and will return whether the model
// converged, and how many iterations it took to converge. If converged == false, the results in mm are meaningless.
//
// To reduce the number of allocations required for repeated use of RunMixtureModel, the input mixture model 'mm' can be reused between calls
// with no modifications necessary.
func RunMixtureModel1D(data []float64, k int, maxIterations int, maxResets int, mm *MixtureModel) (converged bool, iterationsRun int) {
	if len(data) == 0 {
		return
	}
	initMixtureModel(data, k, maxIterations, mm)
	var resets int
	var prevLogLikelihood float64

	for iterationsRun = 0; resets < maxResets && !converged; iterationsRun++ {

		// E step
		prevLogLikelihood = mm.LogLikelihood
		expectation(mm)
		if iterationsRun > 2 && math.Abs(mm.LogLikelihood-prevLogLikelihood) < logProbEpsilon {
			converged = true
		}

		// M step
		maximization(mm)

		// check to see if variance or weights are going to zero which is undesirable as it indicates that
		// the model is stuck overfitting or trying to fit the data to fewer gaussians than intended.
		for i := 0; i < mm.K; i++ {
			switch {
			case mm.Stdev[i] < 1e-04: // most likely caused by overfit of gaussian to a single data point
				fallthrough
			case mm.Weights[i] < 1e-02: // they asked for 2 gaussians, reset and see if you can get better answer
				resets++
				initMixtureModel(data, k, maxIterations, mm)
				iterationsRun = 0
				prevLogLikelihood = 0
				converged = false
				break
			}
		}

		// if max iterations reached, reset with new values and try again
		if iterationsRun == mm.MaxIter {
			resets++
			initMixtureModel(data, k, maxIterations, mm)
			iterationsRun = 0
			prevLogLikelihood = 0
			converged = false
		}
	}
	return
}

// initMixtureModel allocates and fills the MixtureModel struct. Avoids allocation if enough space is already present in input struct.
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
		mm.Means[i] = data[RandIntInRange(0, len(data)-1)]
		mm.Stdev[i] = 1
	}

	if cap(mm.residuals) >= k {
		mm.residuals = mm.residuals[0:k]
		mm.posteriors = mm.posteriors[0:k]
	} else {
		mm.residuals = make([][]float64, k)
		mm.posteriors = make([][]float64, k)
	}

	for i := range mm.residuals {
		if cap(mm.residuals[i]) >= len(data) {
			mm.residuals[i] = mm.residuals[i][0:len(data)]
			mm.posteriors[i] = mm.posteriors[i][0:len(data)]
		} else {
			mm.residuals[i] = make([]float64, len(data))
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

// resetResSum zeros the posteriorSum slice
func resetResSum(mm *MixtureModel) {
	for i := range mm.posteriorsSum {
		mm.posteriorsSum[i] = 0
	}
}

// expectation is the first half of the EM algorithm and determines how well the observed data fit the current model
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
			r = r * r
			mm.residuals[j][i] = r
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

// maximization is the second half of the EM algorithm and generates a new model based on the performance of the previous model
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

	var std, mu float64
	for j := 0; j < mm.K; j++ {
		mu = 0
		std = 0
		for i := range mm.Data {
			mu += mm.posteriors[j][i] * mm.Data[i]
		}

		if mm.posteriorsSum[j] > 0 {
			mm.Means[j] = mu / mm.posteriorsSum[j]
		}

		for i := range mm.Data {
			std += mm.posteriors[j][i] * mm.residuals[j][i]
		}

		if mm.posteriorsSum[j] > 0 {
			std = std / mm.posteriorsSum[j]
		}

		mm.Stdev[j] = math.Sqrt(std)
	}
}
