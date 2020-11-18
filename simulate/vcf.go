package simulate

import (
	"log"
	"math"
	"github.com/vertgenlab/gonomics/numbers"
	"github.com/vertgenlab/gonomics/popgen"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/vcf"
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
		genotype, temp = popgen.SimulateGenotype(alpha, n)
		counts = append(counts, temp)
		//most fields are hardcoded but can be filled in later
		current = &vcf.Vcf{Chr: "chr1", Pos: int64(i+1), Id: ".", Ref: "A", Alt: "T", Qual: 100, Filter: ".", Info: ".", Format: ".", Notes: vcf.GenotypeToStringNew(genotype)}
		vcf.WriteVcf(out, current)
	}

	//log.Printf("testCalc:%f, %f\n", math.Exp(AFSSampleDensity(100, 50, -4, 0, 1, 1e-12)))
	integrateMe := FIntegralComponent(500, 1, -4)
	integrateMeNonLog := FIntegralComponentNonLog(500, 1, -4)
	integrateMeCareful := FIntegralComponentCareful(500, 1, -4)
	integrateMeCarefulLog := FIntegralComponentCarefulLog(500, 1, -4)
	log.Printf("eval:%f, %f, %f, %f\n", math.Exp(integrateMe(0.4)), integrateMeNonLog(0.4), integrateMeCareful(0.4), math.Exp(integrateMeCarefulLog(0.4)))
	log.Printf("integral:%e\n", math.Exp(numbers.AdaptiveSimpsonsLog(integrateMe, 0, 1, 1e-10, 1000)))
	log.Printf("integral:%e\n", numbers.AdaptiveSimpsons(integrateMeNonLog, 0, 1, 1e-10, 1000))
	log.Printf("integral:%e\n", numbers.AdaptiveSimpsons(integrateMeCareful, 0, 1, 1e-10, 1000))
	log.Printf("integral:%e\n", math.Exp(numbers.AdaptiveSimpsonsLog(integrateMeCarefulLog, 0, 1, 1e-10, 1000)))

	var testAlpha, currCalc, ratio float64
	var shouldBeBest float64 = AFSLikelihood(counts, n, alpha)
	for testAlpha = alpha-2; testAlpha <= alpha+3; testAlpha+=0.2 {
		if alpha == 0 {
			currCalc = AFSLikelihood(counts, n, 0.03)
			ratio = math.Exp(currCalc-shouldBeBest)
			log.Printf("like at %f:%e ratio:%f\n", 0.03, AFSLikelihood(counts, n, 0.03))
		} else {
			currCalc = AFSLikelihood(counts, n, testAlpha)
			ratio = math.Exp(currCalc-shouldBeBest)
			log.Printf("like at %f:%e ratio:%f\n", testAlpha, AFSLikelihood(counts, n, testAlpha), ratio)
		}
	}
}

func FIntegralComponentCarefulLog(n int, k int, alpha float64) func(float64) float64 {
        var binomCoeff float64 = numbers.BinomCoefficientLog(n, k)
        return func(p float64) float64 {
                if p == 0 || p == 1 {
                        return math.Inf(-1)
                }
                var ans float64 = binomCoeff + float64(k-1) * math.Log(p) + float64(n-k-1) * math.Log(1.0-p) + math.Log((1.0-math.Exp(-1.0*alpha*(1.0-p))) * 2.0 / (1.0-math.Exp(-1.0*alpha)))
                return ans
        }
}

func FIntegralComponentCareful(n int, k int, alpha float64) func(float64) float64 {
	var binomCoeff float64 = numbers.BinomCoefficientLog(n, k)
        return func(p float64) float64 {
		if p == 0 || p == 1 {
			return 0
		}
                var ans float64 = binomCoeff + float64(k-1) * math.Log(p) + float64(n-k-1) * math.Log(1.0-p) + math.Log((1.0-math.Exp(-1.0*alpha*(1.0-p))) * 2.0 / (1.0-math.Exp(-1.0*alpha)))
		return math.Exp(ans)
        }
}

func FIntegralComponentNonLog(n int, k int, alpha float64) func(float64) float64 {
        return func(p float64) float64 {
		if p == 0 || p == 1 {
                        return 0
                }
		ans := math.Exp(numbers.BinomCoefficientLog(n, k)) * math.Pow(p, float64(k-1)) * math.Pow(1-p, float64(n-k-1)) * (1.0-math.Exp(-1.0*alpha*(1.0-p))) * 2.0 / (1.0-math.Exp(-1.0*alpha))
                return ans
        }
}

func FIntegralComponent(n int, k int, alpha float64) func(float64) float64 {
	var binomCoeff float64 = numbers.BinomCoefficientLog(n, k)
        return func(p float64) float64 {
		if p == 0 || p == 1 {
                        return math.Inf(-1)
                }
                expression := numbers.BinomialExpressionLog(n-2, k-1, p)
		//expressionTwo := float64(k-1) * math.Log(p) + float64(n-k-1) * math.Log(1.0-p)
		//if math.Abs(expression-expressionTwo) > 1e-6 {
		//	log.Fatalf("mismatch expression: n=%d k=%d p=%e %e %e\n", n, k, p, expression, expressionTwo)
		//}
                logPart := math.Log((1.0-math.Exp(-1.0*alpha*(1.0-p))) * 2.0 / (1.0-math.Exp(-1.0*alpha)))
                return numbers.MultiplyLog(binomCoeff,numbers.MultiplyLog(expression, logPart))
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
        constantComponent := 0.0
	//constantComponent := numbers.BinomCoefficientLog(n, k)
        return numbers.MultiplyLog(constantComponent, numbers.AdaptiveSimpsonsLog(f, start, end, accuracy, 1000))
}

//AlleleFrequencyProbability returns the probability of observing i out of n alleles from a stationarity distribution with selection parameter alpha.
func AlleleFrequencyProbability(i int, n int, alpha float64, start, end, accuracy float64) float64 {
        var denom, denomCareful, denomCarefulLog float64 = math.Inf(-1), math.Inf(-1), math.Inf(-1)//denominator begins at -Inf when in log space
	var numer, numerCareful, numerCarefulLog float64
	var ans, ansCareful, ansCarefulLog float64
        //check if n has already been seen
        for j := 1; j < n; j++ {
                denom = numbers.AddLog(denom, AFSSampleDensity(n, j, alpha, start, end, accuracy))
		denomCareful = numbers.AddLog(denomCareful, AFSSampleDensityCareful(n, j, alpha, start, end, accuracy))
		denomCarefulLog = numbers.AddLog(denomCarefulLog, AFSSampleDensityCarefulLog(n, j, alpha, start, end, accuracy))
        }
	numer = AFSSampleDensity(n, i, alpha, start, end, accuracy) // original answer log, integral log
	numerCareful = AFSSampleDensityCareful(n, i, alpha, start, end, accuracy) // answer log, integral normal
	numerCarefulLog = AFSSampleDensityCarefulLog(n, i, alpha, start, end, accuracy) // answer log, integral log
	ans = numbers.DivideLog(numer, denom)
	ansCareful = numbers.DivideLog(numerCareful, denomCareful)
	ansCarefulLog = numbers.DivideLog(numerCarefulLog, denomCarefulLog)
	if math.Abs(ans-ansCareful) > 1e-5 || math.Abs(ans-ansCarefulLog) > 1e-5 { // this is for normal space, but need log space
		log.Printf("%e %e %e\n", ans, ansCareful, ansCarefulLog)
		log.Printf("%d %d %f\n", n, i, alpha)
	}
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
