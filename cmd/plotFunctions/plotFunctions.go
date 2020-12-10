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

func plotFunctions(function string, functionArgs string, left float64, right float64, bins int, outFile string) {
	if function == "AfsStationarity" {
		words := strings.Split(functionArgs, ",")
		if len(words) != 1 {
			log.Fatalf("A stationarity distribution is defined by one parameter, received %d.", len(words))
		}
		alpha := common.StringToFloat64(words[0])
		f := popgen.AFSStationarityClosure(alpha)
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
		fmt.Printf("Unrecognized function: %s.\n", function)
	}
}

func usage() {
	fmt.Print(
		"plotFunctions-returns a tab separated list of function evaluations for plotting functions.\n" +
			"To specify the function, use a function keyword. Then provide the function arguments as a comma separated list in the second argument.\n" +
			"Usage:\n" +
			" plotFunctions functionKeyWord functionArgs leftBound rightBound steps outFile\n" +
			"For AfsStationarity: defined by one parameter (called alpha). ex usage: plotFunctions AfsStationarity 0.001 0.001 0.999 1000 afsPlot.txt\n" + 
			"For Beta: defined by two parameters (called alpha and beta). ex usage: plotFunctions Beta 0.5,0.5 0.001 0.999 1000 betaPlot.txt\n" +
			"For Gamma: defined by two parameters (called alpha and beta). ex usage: plotFunctions Gamma 0.5,0.5 0.001 0.999 1000 gammaPlot.txt\n" +
			"For Normal: defined by two parameters (called mu and sigma). ex usage: plotFunctions Normal 0,0.5 -1 1 1000 normalPlot.txt\n" +
			"options:\n")
	flag.PrintDefaults()
}

func main() {
	var expectedNumArgs int = 6
	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime | log.Lshortfile)
	flag.Parse()

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

	plotFunctions(function, functionArgs, left, right, bins, outFile)
}
