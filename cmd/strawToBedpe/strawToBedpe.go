// Command Group: "Data Conversion"

// Convert HiC contact maps in straw format to bedpe contact peak calls
package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/hic"
	"github.com/vertgenlab/gonomics/numbers"
	"github.com/vertgenlab/gonomics/numbers/fit"
	"log"
	"math"
	"strings"
)

func strawToBedpe(fileList string, outFile string, binSize int) {
	var currStrawChan <-chan hic.Straw
	var words []string
	var binDistance int
	var currDistance float64
	var contactScoreCache [][]int = make([][]int, 2)
	contactScoreCache[0] = make([]int, 1)
	contactScoreCache[1] = make([]int, 1)
	var tmpContactScoreCache [][]int
	var tmpColumn []int

	lines := fileio.Read(fileList)

	for i := range lines {
		words = strings.Split(lines[i], "\t")
		currStrawChan = hic.GoReadToChan(words[0])
		for s := range currStrawChan {

			currDistance = math.Abs(float64(s.Bin1Start) - float64(s.Bin2Start))
			if int(currDistance)%binSize != 0 {
				log.Fatalf("Error: Distance between two straw ends: %v is not a multiple of the bin size: %v.\n", currDistance, binSize)
			}
			binDistance = int(currDistance) / binSize
			// if we need more room in the cache rows, extend the cache
			if binDistance > len(contactScoreCache)-1 {
				tmpContactScoreCache = make([][]int, binDistance+1)
				copy(tmpContactScoreCache, contactScoreCache)
				contactScoreCache = tmpContactScoreCache
			}

			if contactScoreCache[binDistance] == nil {
				contactScoreCache[binDistance] = make([]int, 1)
			}
			// if we need more room in the cache columns, extend the cache
			if s.ContactScore > len(contactScoreCache[binDistance])-1 {
				tmpColumn = make([]int, s.ContactScore+1)
				copy(tmpColumn, contactScoreCache[binDistance])
				contactScoreCache[binDistance] = tmpColumn
			}
			contactScoreCache[binDistance][s.ContactScore]++
		}
	}
	//fmt.Print(contactScoreCache)

	distributions := fitNegativeBinomialsToCountData(contactScoreCache)
	for i := range distributions {
		fmt.Printf("binDistance: %v. R: %v. P: %v.\n", i, distributions[i].R, distributions[i].P)
	}

	var err error
	for i := 20; i < 40; i += 2 {
		out := fileio.EasyCreate(fmt.Sprintf("distFiles/%v.txt", i))
		_, err = fmt.Fprintf(out, "Score\tFrequency\n")
		exception.PanicOnErr(err)
		for j := range contactScoreCache[i] {
			_, err = fmt.Fprintf(out, "%v\t%v\n", j, contactScoreCache[i][j])
			exception.PanicOnErr(err)
		}
		err = out.Close()
		exception.PanicOnErr(err)
	}
}

// NegativeBinomialFit stores the R and P parameter values to specify
// a negative binomial distribution.
type NegativeBinomialFit struct {
	R float64
	P float64
}

// fitNegativeBinomialsToCountData takes in the contactScoreCache, imputes missing
// values at f(0), and fits a negative binomial distribution for each bin distance.
func fitNegativeBinomialsToCountData(cache [][]int) []NegativeBinomialFit {
	var answer = make([]NegativeBinomialFit, len(cache))
	var i int
	var currSumNonZero int
	var currMean, currVariance float64
	var currCount int
	var iteration int
	var couldNotFit bool
	var currDensity float64
	var lastZeroDensity float64 //last iteration's estimate of the density at f(0)

	for currBinDistance := range cache {
		if cache[currBinDistance] == nil {
			cache[currBinDistance] = make([]int, 1)
		}

		/*
			// we initialize contactScoreCache[currBinDistance][0] to the LagrangeInterpolation prediction
			// if we only have two datapoints (f(1) and f(2)), we estimate f(0) by linear interpolation
			// note that if len(contactScoreCache[curBinDistance]) < 3, we cannot interpolate f(0) and cannot call bedpe peaks for these distances.
			if len(cache[currBinDistance]) == 3 {
				cache[currBinDistance][0] = int(fit.LagrangeInterpolation(0.0, [][]float64{
					{1.0, float64(cache[currBinDistance][1])},
					{2.0, float64(cache[currBinDistance][2])},
				},
				))
			} else if len(cache[currBinDistance]) > 3 {
				cache[currBinDistance][0] = int(fit.LagrangeInterpolation(0.0, [][]float64{
					{1.0, float64(cache[currBinDistance][1])},
					{2.0, float64(cache[currBinDistance][2])},
					{3.0, float64(cache[currBinDistance][3])},
				},
				))
			}
		*/

		//if len(cache[currBinDistance]) > 1 {
		//	cache[currBinDistance][0] = cache[currBinDistance][1]
		//}

		// now we fine-tune contactScoreCache[currBinDistance][0] with iterative negative binomial fitting
		currSumNonZero = 0
		currCount = 0

		for i = range cache[currBinDistance] {
			currSumNonZero += cache[currBinDistance][i]
			currCount += cache[currBinDistance][i] * i
		}
		currMean = float64(currCount) / float64(currSumNonZero+cache[currBinDistance][0])

		currVariance = 0
		for i = range cache[currBinDistance] {
			currVariance += float64(cache[currBinDistance][i]) * math.Pow(float64(i)-currMean, 2)
		}
		currVariance = currVariance / (float64(currSumNonZero+cache[currBinDistance][0]) - 1)

		answer[currBinDistance] = NegativeBinomialFit{-1, -1} // initialize to -1 to denote missing data.
		// if we have at least 1000 non-zero entries, we can try to fit  negative binomial to the data
		lastZeroDensity = 0
		currDensity = 1
		iteration = 0
		if currSumNonZero > 1000 {
			for math.Abs(currDensity-lastZeroDensity)/currDensity > 1e-6 && iteration < 100 {
				currCount = 0
				for i = range cache[currBinDistance] {
					currCount += cache[currBinDistance][i] * i
				}
				currMean = float64(currCount) / float64(currSumNonZero+cache[currBinDistance][0])

				currVariance = 0
				for i = range cache[currBinDistance] {
					currVariance += float64(cache[currBinDistance][i]) * math.Pow(float64(i)-currMean, 2)
				}
				currVariance = currVariance / (float64(currSumNonZero+cache[currBinDistance][0]) - 1)
				answer[currBinDistance].R, answer[currBinDistance].P, couldNotFit = fit.NegativeBinomialFromSumStats(currMean, currVariance)
				if couldNotFit {
					answer[currBinDistance].R = -1
					answer[currBinDistance].P = -1
					break
				}
				//fmt.Printf("CurrDistance: %v. Iteration: %v. CurrDensity: %v. LastZeroDensity: %v.\n", currBinDistance, iteration, currDensity, lastZeroDensity)
				lastZeroDensity = currDensity
				currDensity = numbers.NegativeBinomialDist(0, answer[currBinDistance].R, answer[currBinDistance].P)
				cache[currBinDistance][0] = int(currDensity * float64(currSumNonZero+cache[currBinDistance][0]))
				iteration++
			}
		}
	}

	return answer
}

func usage() {
	fmt.Print(
		"strawToBedpe - Convert HiC contact maps in straw format to bedpe contact peak calls.\n" +
			"The user must provide a txt file with the following tab separated fields:\n" +
			"file1.straw\tchrom\n" +
			"In the line above, the first field specifies a straw file, and the second field specifies\n" +
			"the chromosome name.\n" +
			"Inter-chromosomal straw files (where aChrom != bChrom), are not currently supported.\n" +
			"\n" +
			"Usage:\n" +
			"strawToBedpe [options] fileListAndChromNames.txt out.bedpe\n" +
			"options:\n")
	flag.PrintDefaults()
}

func main() {
	var expectedNumArgs int = 2
	var binSize *int = flag.Int("binSize", 5000, "Specify the binSize of input straw data.")

	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime | log.Lshortfile)
	flag.Parse()

	if len(flag.Args()) != expectedNumArgs {
		flag.Usage()
		log.Fatalf("Error: expecting %d arguments, but got %d\n",
			expectedNumArgs, len(flag.Args()))
	}

	fileListAndChromFile := flag.Arg(0)
	outFile := flag.Arg(1)

	strawToBedpe(fileListAndChromFile, outFile, *binSize)
}
