package simulate

import (
	"fmt"
	"github.com/vertgenlab/gonomics/common"
	"github.com/vertgenlab/gonomics/fileio"
	"math"
	"strings"
	"testing"
)

/*
func getStation(n int, k int, alpha float64) func(float64) float64 {
        var binomCoeff float64 = numbers.BinomCoefficientLog(n, k)
        return func(p float64) float64 {
                expression := numbers.BinomialExpressionLog(n-2, k-1, p)
                logPart := math.Log((1.0 - math.Exp(-1.0*alpha*(1.0-p))) * 2.0 / (1.0 - math.Exp(-1.0*alpha)))
                return numbers.MultiplyLog(binomCoeff, numbers.MultiplyLog(expression, logPart))
        }
}

func binomInIntegral(n int, k int, alpha float64, start float64, end float64, accuracy float64) float64 {
        f := getStation(n, k, alpha)
        integralResult := numbers.AdaptiveSimpsonsLog(f, start, end, accuracy, 1000)
        return integralResult
}

func binomInIntegralParts(n int, k int, alpha float64, start float64, end float64, accuracy float64) float64 {
	var switchPoint float64 = float64(k)/float64(n)
	if switchPoint <= start || switchPoint >= end {
		log.Fatal("Error: confused about switch point\n")
	}
	f := getStation(n, k, alpha)
        partOne := numbers.AdaptiveSimpsonsLog(f, start, switchPoint, accuracy, 1000)
	partTwo := numbers.AdaptiveSimpsonsLog(f, switchPoint, end, accuracy, 1000)
        return numbers.AddLog(partOne, partTwo)
}
*/

func TestIntegrateOverAfs(t *testing.T) {
	f := fileio.EasyOpen("testdata/testVal.tsv")
	defer f.Close()

	var line string
	var done bool
	var n, k, good, bad int
	var alpha, expected, checkMe, relativeError float64
	var words []string

	for line, done = fileio.EasyNextRealLine(f); !done; line, done = fileio.EasyNextRealLine(f) {
		words = strings.Split(line, "\t")
		n = common.StringToInt(words[0])
		k = common.StringToInt(words[1])
		alpha = common.StringToFloat64(words[2])
		expected = common.StringToFloat64(words[3])

		//checkMe = math.Exp(binomInIntegral(n, k, alpha, 0, 1, 1e-10))
		checkMe = math.Exp(binomInIntegralParts(n, k, alpha, 1e-10))

		relativeError = math.Abs(expected-checkMe) / expected
		if relativeError > 0.5 {
			fmt.Printf("bad: %d %d %f %f (%e) %f (%e) %f\n", n, k, alpha, checkMe, checkMe, expected, expected, relativeError)
			bad++
			/*if bad > 10 {
				log.Fatal("too many misses\n")
			}*/
		} else {
			good++
		}
	}
	fmt.Printf("Good:%d Bad:%d\n", good, bad)
}
