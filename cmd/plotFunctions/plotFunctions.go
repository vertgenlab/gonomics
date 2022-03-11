// Command Group: "Statistics & Population Genetics"

package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/common"
	"github.com/vertgenlab/gonomics/numbers"
	"github.com/vertgenlab/gonomics/popgen"
	"log"
	"strings"
)

func plotContinuousFunctions(function string, functionArgs string, left float64, right float64, bins int, outFile string) {
	if function == "AfsStationarity" {
		words := strings.Split(functionArgs, ",")
		if len(words) != 1 {
			log.Fatalf("A stationarity distribution is defined by one parameter, received %d.", len(words))
		}
		alpha := common.StringToFloat64(words[0])
		f := popgen.AfsStationarityClosure(alpha)
		numbers.Plot(f, left, right, bins, outFile)
	} else if function == "Beta" {
		words := strings.Split(functionArgs, ",")
		if len(words) != 2 {
			log.Fatalf("A beta distribution is defined by two parameters, received %d.", len(words))
		}
		alpha := common.StringToFloat64(words[0])
		beta := common.StringToFloat64(words[1])
		f := numbers.BetaClosure(alpha, beta)
		numbers.Plot(f, left, right, bins, outFile)
	} else if function == "Gamma" {
		words := strings.Split(functionArgs, ",")
		if len(words) != 2 {
			log.Fatalf("A gamma distribution is defined by two parameters, received %d.", len(words))
		}
		alpha := common.StringToFloat64(words[0])
		beta := common.StringToFloat64(words[1])
		f := numbers.GammaClosure(alpha, beta)
		numbers.Plot(f, left, right, bins, outFile)
	} else if function == "Normal" {
		words := strings.Split(functionArgs, ",")
		if len(words) != 2 {
			log.Fatalf("a normal distribution is defined by two parameters, received %d.", len(words))
		}
		mu := common.StringToFloat64(words[0])
		sigma := common.StringToFloat64(words[1])
		f := numbers.NormalClosure(mu, sigma)
		numbers.Plot(f, left, right, bins, outFile)
	} else { //here you can add more else ifs to add additional functions for plotting
		log.Fatalf("Unrecognized function: %s.\n", function)
	}
}

func usage() {
	fmt.Print(
		"plotFunctions-returns a tab separated list of function evaluations for plotting functions. f(x) is printed in logSpace unless stated otherwise.\n" +
			"To specify the function, use a function keyword. Then provide the function arguments as a comma separated list in the second argument.\n" +
			"Usage:\n" +
			" plotFunctions functionKeyWord functionArgs leftBound rightBound steps outFile\n" +
			"For AfsStationarity: defined by one parameter (called alpha). ex usage: plotFunctions AfsStationarity 0.001 0.001 0.999 1000 afsPlot.txt\n" +
			"For AfsF: defined by three parameters (alpha, n, and integralError: where alpha is the selection parameter and n is the total number of individuals. ex. Usage: AfsF 0.5,200,1e-7 afsFPlot.txt\n" +
			"For AfsProbabilityAncestral: defined by two parameters (alpha and n: where alpha is the selection parameter and n is the total number of individuals. ex. Usage: AfsProbabilityAncestral 0.5,200 afsProbAncestralPlot.txt\n" +
			"For AfsProbabilityDerived: defined by two parameters (alpha and n: where alpha is the selection parameter and n is the total number of individuals. ex. Usage: AfsProbabilityDerived 0.5,200 afsProbDerivedPlot.txt\n" +
			"For AfsProbability: defined by two parameters (alpha and n: where alpha is the selection parameter and n is the total number of individuals. ex. Usage: AfsProbability 0.5,200,1e-7 afsProbabilityPlot.txt\n" +
			"For AscertainmentProbabilityDerived: defined by two parameters (n and d), where n is the total number of individuals and d is the ascertainment subset size. ex. Usage: AscertainmentProbabilityDerived 1002,1 ascertainmentProbabilityDerivedPlot.txt\n" +
			"For AscertainmentProbabilityAncestral: defined by two parameters (n and d), where n is the total number of individuals and d is the ascertainment subset size. ex. Usage: AscertainmentProbabilityAncestral 1002,1 ascertainmentProbabilityAncestralPlot.txt\n" +
			"For AncestralAscertainmentDenominator: defined by four parameters (n, number of total alleles; d, ascertainment subset size; alpha, selection coefficient; integralError). ex.Usage: AncestralAscertainmentDenominator 1002,1,0.03,1e-5 ancestralAscertainmentDenominatorPlot.txt\n" +
			"For DerivedAscertainmentDenominator: defined by four parameters (n, number of total alleles; d, ascertainment subset size; alpha, selection coefficient; integralError). ex.Usage: DerivedAscertainmentDenominator 1002,1,0.03,1e-5 derivedAscertainmentDenominatorPlot.txt \n" +
			"For Beta: defined by two parameters (called alpha and beta). ex usage: plotFunctions Beta 0.5,0.5 0.001 0.999 1000 betaPlot.txt\n" +
			"For Gamma: defined by two parameters (called alpha and beta). ex usage: plotFunctions Gamma 0.5,0.5 0.001 0.999 1000 gammaPlot.txt\n" +
			"For Normal: defined by two parameters (called mu and sigma). ex usage: plotFunctions Normal 0,0.5 -1 1 1000 normalPlot.txt\n" +
			"For ChooseN: defined by 1 parameter (called n). Prints the binomial coefficient. ex usage: plotFunctions ChooseN 40 chooseNPlot.txt\n" +
			"options:\n")
	flag.PrintDefaults()
}

func main() {
	var expectedNumArgs int
	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime | log.Lshortfile)
	flag.Parse()

	if len(flag.Args()) == 0 {
		flag.Usage()
		log.Fatalf("No arguments inputed, help message printed.\n")
	}

	if flag.Arg(0) == "AfsProbability" {

		expectedNumArgs = 3
		if len(flag.Args()) != expectedNumArgs {
			flag.Usage()
			log.Fatalf("Error: expecting %d arguments, but got %d\n", expectedNumArgs, len(flag.Args()))
		}
		functionArgs := flag.Arg(1)
		outFile := flag.Arg(2)
		words := strings.Split(functionArgs, ",")
		if len(words) != 3 {
			log.Fatalf("An allele frequency probability mass function is defined by three parameters, received %d.", len(words))
		}
		alpha := common.StringToFloat64(words[0])
		n := common.StringToInt(words[1])
		integralError := common.StringToFloat64(words[2])
		popgen.PlotAfsPmf(alpha, n, outFile, integralError, false, false)
	} else if flag.Arg(0) == "AfsProbabilityAncestral" {
		expectedNumArgs = 3
		if len(flag.Args()) != 3 {
			flag.Usage()
			log.Fatalf("Error: expecting %d arguments, but got %d\n", expectedNumArgs, len(flag.Args()))
		}
		functionArgs := flag.Arg(1)
		outFile := flag.Arg(2)
		words := strings.Split(functionArgs, ",")
		if len(words) != expectedNumArgs {
			log.Fatalf("An allele frequency probability mass function is defined by three parameters, received %d.", len(words))
		}
		alpha := common.StringToFloat64(words[0])
		n := common.StringToInt(words[1])
		integralError := common.StringToFloat64(words[2])
		popgen.PlotAfsPmf(alpha, n, outFile, integralError, false, true)
	} else if flag.Arg(0) == "AfsProbabilityDerived" {
		expectedNumArgs = 3
		if len(flag.Args()) != expectedNumArgs {
			flag.Usage()
			log.Fatalf("Error: expecting %d arguments, but got %d\n", expectedNumArgs, len(flag.Args()))
		}
		functionArgs := flag.Arg(1)
		outFile := flag.Arg(2)
		words := strings.Split(functionArgs, ",")
		if len(words) != 3 {
			log.Fatalf("An allele frequency probability mass function is defined by three parameters, received %d.", len(words))
		}
		alpha := common.StringToFloat64(words[0])
		n := common.StringToInt(words[1])
		integralError := common.StringToFloat64(words[2])
		popgen.PlotAfsPmf(alpha, n, outFile, integralError, true, false)
	} else if flag.Arg(0) == "AscertainmentProbabilityDerived" {
		expectedNumArgs = 3
		if len(flag.Args()) != expectedNumArgs {
			flag.Usage()
			log.Fatalf("Error: expecting %d arguments, but got %d\n", expectedNumArgs, len(flag.Args()))
		}
		functionArgs := flag.Arg(1)
		outFile := flag.Arg(2)
		words := strings.Split(functionArgs, ",")
		if len(words) != 2 {
			log.Fatalf("An ascertainment probability function is defined by two parameters, received %d.", len(words))
		}
		n := common.StringToInt(words[0])
		d := common.StringToInt(words[1])
		popgen.PlotDerivedAscertainmentProbability(outFile, n, d)
	} else if flag.Arg(0) == "AscertainmentProbabilityAncestral" {
		expectedNumArgs = 3
		if len(flag.Args()) != expectedNumArgs {
			flag.Usage()
			log.Fatalf("Error: expecting %d arguments, but got %d\n", expectedNumArgs, len(flag.Args()))
		}
		functionArgs := flag.Arg(1)
		outFile := flag.Arg(2)
		words := strings.Split(functionArgs, ",")
		if len(words) != 2 {
			log.Fatalf("An allele frequency probability mass function is defined by two parameters, received %d.", len(words))
		}
		n := common.StringToInt(words[0])
		d := common.StringToInt(words[1])
		popgen.PlotAncestralAscertainmentProbability(outFile, n, d)
	} else if flag.Arg(0) == "AncestralAscertainmentDenominator" {
		expectedNumArgs = 3
		if len(flag.Args()) != expectedNumArgs {
			flag.Usage()
			log.Fatalf("Error: expecting %d arguments, but got %d\n", expectedNumArgs, len(flag.Args()))
		}
		functionArgs := flag.Arg(1)
		outFile := flag.Arg(2)
		words := strings.Split(functionArgs, ",")
		if len(words) != 4 {
			log.Fatalf("An ascertainment denominator is defined by four parameters, received %d.", len(words))
		}
		n := common.StringToInt(words[0])
		d := common.StringToInt(words[1])
		alpha := common.StringToFloat64(words[2])
		integralError := common.StringToFloat64(words[3])
		popgen.PlotAncestralAscertainmentDenominator(outFile, n, d, alpha, integralError)
	} else if flag.Arg(0) == "DerivedAscertainmentDenominator" {
		expectedNumArgs = 3
		if len(flag.Args()) != expectedNumArgs {
			flag.Usage()
			log.Fatalf("Error: expecting %d arguments, but got %d\n", expectedNumArgs, len(flag.Args()))
		}
		functionArgs := flag.Arg(1)
		outFile := flag.Arg(2)
		words := strings.Split(functionArgs, ",")
		if len(words) != 4 {
			log.Fatalf("An ascertainment denominator is defined by four parameters, received %d.", len(words))
		}
		n := common.StringToInt(words[0])
		d := common.StringToInt(words[1])
		alpha := common.StringToFloat64(words[2])
		integralError := common.StringToFloat64(words[3])
		popgen.PlotDerivedAscertainmentDenominator(outFile, n, d, alpha, integralError)
	} else if flag.Arg(0) == "ChooseN" {
		expectedNumArgs = 3
		if len(flag.Args()) != expectedNumArgs {
			flag.Usage()
			log.Fatalf("Error: expecting %d arguments, but got %d\n", expectedNumArgs, len(flag.Args()))
		}
		functionArgs := flag.Arg(1)
		outFile := flag.Arg(2)
		n := common.StringToInt(functionArgs)
		numbers.PlotBinomCoefficient(n, outFile)
	} else if flag.Arg(0) == "AfsF" {
		expectedNumArgs = 3
		if len(flag.Args()) != expectedNumArgs {
			flag.Usage()
			log.Fatalf("Error: expecting %d arguments, but got %d\n", expectedNumArgs, len(flag.Args()))
		}
		functionArgs := flag.Arg(1)
		outFile := flag.Arg(2)
		words := strings.Split(functionArgs, ",")
		if len(words) != 3 {
			log.Fatalf("An allele frequency F function is defined by three parameters, received %d.", len(words))
		}
		alpha := common.StringToFloat64(words[0])
		n := common.StringToInt(words[1])
		error := common.StringToFloat64(words[2])
		popgen.PlotAfsF(alpha, n, outFile, error)
	} else {
		expectedNumArgs = 6
		if len(flag.Args()) != expectedNumArgs {
			flag.Usage()
			log.Fatalf("Error: expecting %d arguments, but got %d\n",
				expectedNumArgs, len(flag.Args()))
		}

		function := flag.Arg(0)
		functionArgs := flag.Arg(1)
		left := common.StringToFloat64(flag.Arg(2))
		right := common.StringToFloat64(flag.Arg(3))
		bins := common.StringToInt(flag.Arg(4))
		outFile := flag.Arg(5)

		plotContinuousFunctions(function, functionArgs, left, right, bins, outFile)
	}
}
