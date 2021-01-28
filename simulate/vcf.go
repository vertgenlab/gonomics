package simulate

import (
	"github.com/vertgenlab/gonomics/fileio"
	//"github.com/vertgenlab/gonomics/numbers"
	"github.com/vertgenlab/gonomics/popgen"
	"github.com/vertgenlab/gonomics/vcf"
	"log"
	"math"
)

//SimulateVCF generates simulated VCF data. This currently is a skeleton that can be expanded on with additional functions. Currently supports generating gVCF annotations from
//SimulateAFS in the popgen package. Random locations can be generated with simulateBed, and random mutations can be picked with Christi's simulate code.
//maybe we can combine these?
func SimulateVcf(alpha float64, n int, k int, outFile string) {
	out := fileio.EasyCreate(outFile)
	defer out.Close()
	var genotype []vcf.GenomeSample
	var counts []int
	var temp int

	var current *vcf.Vcf
	//for each segregating site, we make a vcf entry and write out
	for i := 0; i < k; i++ {
		//		genotype = popgen.SimulateGenotype(alpha, n)
		genotype, temp = popgen.SimulateGenotype(alpha, n)
		counts = append(counts, temp)
		//most fields are hardcoded but can be filled in later
		current = &vcf.Vcf{Chr: "chr1", Pos: i + 1, Id: ".", Ref: "A", Alt: []string{"T"}, Qual: 100, Filter: ".", Info: ".", Format: []string{"GT"}, Samples: genotype}

		vcf.WriteVcf(out, current)
	}

	/*
	integrateMe := FIntegralComponent(500, 1, -4)
	integrateMeNonLog := FIntegralComponentNonLog(500, 1, -4)
	integrateMeCareful := FIntegralComponentCareful(500, 1, -4)
	integrateMeCarefulLog := FIntegralComponentCarefulLog(500, 1, -4) 
	log.Printf("eval:%e, %e, %e, %e\n", math.Exp(integrateMe(0.4)), integrateMeNonLog(0.4), integrateMeCareful(0.4), math.Exp(integrateMeCarefulLog(0.4)))
	log.Printf("integral:%e\n", math.Exp(numbers.AdaptiveSimpsonsLog(integrateMe, 0, 1, 1e-10, 1000)))
	log.Printf("integral:%e\n", numbers.AdaptiveSimpsons(integrateMeNonLog, 0, 1, 1e-10, 1000))
	log.Printf("integral:%e\n", numbers.AdaptiveSimpsons(integrateMeCareful, 0, 1, 1e-10, 1000))
	log.Printf("integral:%e\n", math.Exp(numbers.AdaptiveSimpsonsLog(integrateMeCarefulLog, 0, 1, 1e-10, 1000)))
	//log.Fatalf("DONE")

	//var outLier float64 = AFSLikelihood(counts, n, 1)
	//var outLierBinomInIntegral float64 = AFSLikelihoodBinomInIntegral(counts, n, 1)
	//log.Printf("like at outlier: %e. With binomInIntegral: %e.\n", outLier, outLierBinomInIntegral)
	var outLierDensity float64 = AFSSampleDensity(100, 3, 1, 0, 1, 1e-10)
	var outLierDensityBinomInIntegral float64 = AFSSampleDensityBinomInIntegral(100, 3, 1, 0, 1, 1e-10)
	log.Printf("density binom outside: %e. density binom inside: %e.\n", outLierDensity, outLierDensityBinomInIntegral)

	log.Fatalf("DONE")
	*/

	var testAlpha, currCalc, ratio float64
	var shouldBeBest float64
	var accuracy = 1e-8

	//shouldBeBest = AFSLikelihood(counts, n, alpha)
	shouldBeBest = AFSLikelihoodNew(counts, n, alpha, accuracy)
	for testAlpha = alpha - 4; testAlpha <= alpha+4; testAlpha += 0.1 {
		if testAlpha == 0 {
			log.Fatal("Error: right now alpha can not be zero\n")
		} else {
			//currCalc = AFSLikelihood(counts, n, testAlpha)
			currCalc = AFSLikelihoodNew(counts, n, testAlpha, accuracy)
			ratio = math.Exp(currCalc - shouldBeBest)
			log.Printf("like at %f:%e ratio:%f\n", testAlpha, currCalc, ratio)
		}
	}
}
/*
func FIntegralComponentCarefulLog(n int, k int, alpha float64) func(float64) float64 {
	var binomCoeff float64 = numbers.BinomCoefficientLog(n, k)
	return func(p float64) float64 {
		//var ans float64 = binomCoeff + float64(k-1)*math.Log(p) + float64(n-k-1)*math.Log(1.0-p) + math.Log((1.0-math.Exp(-1.0*alpha*(1.0-p)))*2.0/(1.0-math.Exp(-1.0*alpha)))
		var ans float64 = binomCoeff + numbers.LogPow(p, float64(k-1)) + numbers.LogPow(1.0 - p, float64(n-k-1)) + math.Log((1.0-math.Exp(-1.0*alpha*(1.0-p)))*2.0/(1.0-math.Exp(-1.0*alpha)))
		return ans
	}
}

func FIntegralComponentCareful(n int, k int, alpha float64) func(float64) float64 {
	var binomCoeff float64 = numbers.BinomCoefficientLog(n, k)
	return func(p float64) float64 {
		//carefulExpression := numbers.LogPow(float64(k-1), p) + numbers.LogPow(float64(n-k-1), 1.0-p)
		carefulExpression := numbers.MultiplyLog(numbers.LogPow(p, float64(k-1)), numbers.LogPow(1.0-p, float64(n-k-1)))
		//functionExpression := numbers.BinomialExpressionLog(n-2, k-1, p)
		//log.Printf("n: %v. k: %v. alpha: %f. p: %v. carefulExpression: %f. functionExpression: %f.\n", n, k, alpha, p, carefulExpression, functionExpression)
		var ans float64 = binomCoeff + carefulExpression + math.Log((1.0-math.Exp(-1.0*alpha*(1.0-p)))*2.0/(1.0-math.Exp(-1.0*alpha)))
		//var ans float64 = binomCoeff + float64(k-1) * math.Log(p) + float64(n-k-1) * math.Log(1.0-p) + math.Log((1.0-math.Exp(-1.0*alpha*(1.0-p)))*2.0/(1.0-math.Exp(-1.0*alpha)))
		return math.Exp(ans)
	}
}

func FIntegralComponentNonLog(n int, k int, alpha float64) func(float64) float64 {
	var binomCoeff float64 = numbers.BinomCoefficientLog(n, k)
	return func(p float64) float64 {
		ans := math.Exp(binomCoeff) * math.Pow(p, float64(k-1)) * math.Pow(1-p, float64(n-k-1)) * (1.0 - math.Exp(-1.0*alpha*(1.0-p))) * 2.0 / (1.0 - math.Exp(-1.0*alpha))
		return ans
	}
}

func FIntegralComponent(n int, k int, alpha float64) func(float64) float64 {
	//var binomCoeff float64 = numbers.BinomCoefficientLog(n, k)
	return func(p float64) float64 {
		expression := numbers.BinomialExpressionLog(n-2, k-1, p)
		//expressionTwo := float64(k-1) * math.Log(p) + float64(n-k-1) * math.Log(1.0-p)
		//if math.Abs(expression-expressionTwo) > 1e-6 {
		//	log.Fatalf("mismatch expression: n=%d k=%d p=%e %e %e\n", n, k, p, expression, expressionTwo)
		//}
		logPart := math.Log((1.0 - math.Exp(-1.0*alpha*(1.0-p))) * 2.0 / (1.0 - math.Exp(-1.0*alpha)))
		//return numbers.MultiplyLog(binomCoeff, numbers.MultiplyLog(expression, logPart))
		return numbers.MultiplyLog(expression, logPart)
		//return numbers.MultiplyLog(expression, logPart)
	}
}

func FIntegralComponentBinomInIntegral(n int, k int, alpha float64) func(float64) float64 {
	return func(p float64) float64 {
		var binomCoeff float64 = numbers.BinomCoefficientLog(n, k)
		expression := numbers.BinomialExpressionLog(n-2, k-1, p)
		//expressionTwo := float64(k-1) * math.Log(p) + float64(n-k-1) * math.Log(1.0-p)
		//if math.Abs(expression-expressionTwo) > 1e-6 {
		//	log.Fatalf("mismatch expression: n=%d k=%d p=%e %e %e\n", n, k, p, expression, expressionTwo)
		//}
		logPart := math.Log((1.0 - math.Exp(-1.0*alpha*(1.0-p))) * 2.0 / (1.0 - math.Exp(-1.0*alpha)))
		return numbers.MultiplyLog(binomCoeff, numbers.MultiplyLog(expression, logPart))
		//return numbers.MultiplyLog(expression, logPart)
		//return numbers.MultiplyLog(expression, logPart)
	}
}

func ficBinomCoeff(n int, k int, alpha float64, binomCoeff float64) func(float64) float64 {
	return func(p float64) float64 {
		expression := numbers.BinomialExpressionLog(n-2, k-1, p)
		logPart := math.Log((1.0 - math.Exp(-1.0*alpha*(1.0-p))) * 2.0 / (1.0 - math.Exp(-1.0*alpha)))
		return numbers.MultiplyLog(binomCoeff, numbers.MultiplyLog(expression, logPart))
		//return numbers.MultiplyLog(expression, logPart)
	}
}

func fic(n int, k int, alpha float64, binomCoeff float64) func(float64) float64 {
	return func(p float64) float64 {
		expression := numbers.BinomialExpressionLog(n-2, k-1, p)
		logPart := math.Log((1.0 - math.Exp(-1.0*alpha*(1.0-p))) * 2.0 / (1.0 - math.Exp(-1.0*alpha)))
		return numbers.MultiplyLog(binomCoeff, numbers.MultiplyLog(expression, logPart))
		//return numbers.MultiplyLog(expression, logPart)
	}
}

func AFSSampleDensityCarefulLog(n int, k int, alpha float64, start float64, end float64, accuracy float64) float64 {
	//DEBUG: log.Printf("n: %d. k: %d. alpha: %v.", n, k, alpha)
	f := FIntegralComponentCarefulLog(n, k, alpha)
	return numbers.AdaptiveSimpsonsLog(f, start, end, accuracy, 1000)
}

func AFSSampleDensityCareful(n int, k int, alpha float64, start float64, end float64, accuracy float64) float64 {
	//DEBUG: log.Printf("n: %d. k: %d. alpha: %v.", n, k, alpha)
	f := FIntegralComponentCareful(n, k, alpha)
	return math.Log(numbers.AdaptiveSimpsons(f, start, end, accuracy, 1000))
}

func AFSSampleDensity(n int, k int, alpha float64, start float64, end float64, accuracy float64) float64 {
	//DEBUG: log.Printf("n: %d. k: %d. alpha: %v.", n, k, alpha)
	f := FIntegralComponent(n, k, alpha)
	//constantComponent := 0.0
	constantComponent := numbers.BinomCoefficientLog(n, k)
	integralResult := numbers.AdaptiveSimpsonsLog(f, start, end, accuracy, 1000)
	log.Printf("AFSSAmpleDensity. integralResult: %e. constantComponent: %e.\n", math.Exp(integralResult), math.Exp(constantComponent))
	return numbers.MultiplyLog(constantComponent, integralResult)
}

func AFSSampleDensityBinomInIntegral(n int, k int, alpha float64, start float64, end float64, accuracy float64) float64 {
	//DEBUG: log.Printf("n: %d. k: %d. alpha: %v.", n, k, alpha)
	f := FIntegralComponentBinomInIntegral(n, k, alpha)
	constantComponent := 0.0
	integralResult := numbers.AdaptiveSimpsonsLog(f, start, end, accuracy, 1000)
	log.Printf("AFSSampleDensityBinomInIntegral. integralResult: %e. constantComponent: %e.\n", math.Exp(integralResult), math.Exp(constantComponent))
	//constantComponent := numbers.BinomCoefficientLog(n, k)
	return numbers.MultiplyLog(constantComponent, integralResult)
}

//AlleleFrequencyProbability returns the probability of observing i out of n alleles from a stationarity distribution with selection parameter alpha.
func AlleleFrequencyProbability(i int, n int, alpha float64, start, end, accuracy float64) float64 {
	var denom float64 = math.Inf(-1)
	var numer float64
	//var denom, denomCareful, denomCarefulLog float64 = math.Inf(-1), math.Inf(-1), math.Inf(-1) //denominator begins at -Inf when in log space
	//var numer, numerCareful, numerCarefulLog float64
	//var ans, ansCareful, ansCarefulLog float64
	//check if n has already been seen
	for j := 1; j < n; j++ {
		denom = numbers.AddLog(denom, AFSSampleDensity(n, j, alpha, start, end, accuracy))
		//denomCareful = numbers.AddLog(denomCareful, AFSSampleDensityCareful(n, j, alpha, start, end, accuracy))
		//denomCarefulLog = numbers.AddLog(denomCarefulLog, AFSSampleDensityCarefulLog(n, j, alpha, start, end, accuracy))
	}
	numer = AFSSampleDensity(n, i, alpha, start, end, accuracy)                     // original answer log, integral log
	//numerCareful = AFSSampleDensityCareful(n, i, alpha, start, end, accuracy)       // answer log, integral normal
	//numerCarefulLog = AFSSampleDensityCarefulLog(n, i, alpha, start, end, accuracy) // answer log, integral log
	log.Printf("For regular: i: %v. numer: %e. denom: %e.\n", i, numer, denom)
	return numbers.DivideLog(numer, denom)
}

func AlleleFrequencyProbabilityBinomInIntegral(i int, n int, alpha float64, start, end, accuracy float64) float64 {
	var denom float64 = math.Inf(-1)
	var numer float64
	//var denom, denomCareful, denomCarefulLog float64 = math.Inf(-1), math.Inf(-1), math.Inf(-1) //denominator begins at -Inf when in log space
	//var numer, numerCareful, numerCarefulLog float64
	//var ans, ansCareful, ansCarefulLog float64
	//check if n has already been seen
	for j := 1; j < n; j++ {
		denom = numbers.AddLog(denom, AFSSampleDensity(n, j, alpha, start, end, accuracy))
		//denomCareful = numbers.AddLog(denomCareful, AFSSampleDensityCareful(n, j, alpha, start, end, accuracy))
		//denomCarefulLog = numbers.AddLog(denomCarefulLog, AFSSampleDensityCarefulLog(n, j, alpha, start, end, accuracy))
	}
	numer = AFSSampleDensityBinomInIntegral(n, i, alpha, start, end, accuracy)                     // original answer log, integral log
	//numerCareful = AFSSampleDensityCareful(n, i, alpha, start, end, accuracy)       // answer log, integral normal
	//numerCarefulLog = AFSSampleDensityCarefulLog(n, i, alpha, start, end, accuracy) // answer log, integral log
	log.Printf("For BinomInIntegral: i: %v. numer: %e. denom: %e.\n", i, numer, denom)
	return numbers.DivideLog(numer, denom)
}

//func AFSLikelihood(count []int, totalAlleles int, alpha float64, start, end, accuracy float64) float64 {
func AFSLikelihood(count []int, totalAlleles int, alpha float64) float64 {
	var answer float64 = 0.0
	for j := 0; j < len(count); j++ {
		answer = numbers.MultiplyLog(answer, AlleleFrequencyProbability(count[j], totalAlleles, alpha, 0, 1, 1e-10))
	}
	return answer
}

func AFSLikelihoodBinomInIntegral(count []int, totalAlleles int, alpha float64) float64 {
	var answer float64 = 0.0
	for j := 0; j < len(count); j++ {
		answer = numbers.MultiplyLog(answer, AlleleFrequencyProbabilityBinomInIntegral(count[j], totalAlleles, alpha, 0, 1, 1e-10))
	}
	return answer
}
*/
