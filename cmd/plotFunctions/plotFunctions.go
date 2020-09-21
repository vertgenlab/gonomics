package main

import (
	"github.com/vertgenlab/gonomics/numbers"
	"github.com/vertgenlab/gonomics/popgen"
	"github.com/vertgenlab/gonomics/common"
	"log"
	"fmt"
	"flag"
	"math"
)

func plotFunctions(function string, left float64, right float64, bins int, outFile string, alpha float64) {
	if function == "AfsStationarity" {
		numbers.Plot(popgen.AFSStationarityClosure(alpha), left, right, bins, outFile)
	} else {//here you can add more else ifs to add additional functions for plotting
		fmt.Printf("Unrecognized function: %s.\n", function)
	}
}

func usage() {
	fmt.Print(
		"plotFunctions-returns a tab separated list of function evaluations for plotting functions.\n" +
		"To specify the function, use a function keyword. Currently only 'AfsStationarity' is supported.\n" + 
			"Usage:\n" +
			" plotFunctions functionKeyWord leftBound rightBound numberOfBins outFile\n" +
			"options:\n")
	flag.PrintDefaults()
}

func main() {
	var expectedNumArgs int = 5
	var alpha *float64 = flag.Float64("alpha", math.Inf(-1), "Specifies the strength of selection for AfsStationarity plots.")
	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime | log.Lshortfile)
	flag.Parse()
	
	if len(flag.Args()) != expectedNumArgs {
		flag.Usage()
		log.Fatalf("Error: expecting %d arguments, but got %d\n",
			expectedNumArgs, len(flag.Args()))
	}

	function := flag.Arg(0)
	left := common.StringToFloat64(flag.Arg(1))
	right := common.StringToFloat64(flag.Arg(2))
	bins := common.StringToInt(flag.Arg(3))
	outFile := flag.Arg(4)

	plotFunctions(function, left, right, bins, outFile, *alpha)
}
