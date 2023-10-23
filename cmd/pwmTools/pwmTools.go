package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/numbers"
	"log"
)

func usage() {
	fmt.Printf(
		"pwmTools - a collection of tools for manipulating position matrices.\n" +
			"Usage:\n" +
			"pwmTools filter in.pwm out.pwm\n" +
			"OR\n" +
			"pwmTools info in.pwm out.txt\n" +
			"options:\n")
	flag.PrintDefaults()
}

func filterUsage() {
	fmt.Printf("pwmTools filter - a tool for filtering PWM records.\n" +
		"Usage:\n" +
		"pwmTools filter in.pwm out.pwm\n" +
		"options:\n")
	flag.PrintDefaults()
}

func infoUsage() {
	fmt.Printf("pwmTools info - a tool for reporting tabular statistics on PWM files.\n" +
		"Usage:\n" +
		"pwmTools info in.pwm out.txt\n" +
		"options:\n")
	flag.PrintDefaults()
}

func main() {
	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime | log.Lshortfile)
	flag.Parse()

	if len(flag.Args()) < 1 {
		flag.Usage()
		log.Fatalf("Error: user must specify a pwmTools subcommand.\n")
	}

	switch flag.Arg(0) {
	case "filter":
		parseFilterArgs()
	case "info":
		parseInfoArgs()
	default:
		flag.Usage()
		log.Fatalf("Error: unrecognized subcommand: %v.\n", flag.Arg(0))
	}
}

func parseFilterArgs() {
	var expectedNumArgs int = 2
	var minLength *int = flag.Int("minLength", 0, "Specify the minimum length of PWMs to be retained in the output.")
	var maxLength *int = flag.Int("maxLength", numbers.MaxInt, "Specify the maximum length of PWMs to be retained in the output.")
	var matrixType *string = flag.String("matrixType", "Frequency", "Specify the type of position matrix. May be one of the following: Frequency, Probability, Weight.")
	flag.Usage = filterUsage
	log.SetFlags(log.Ldate | log.Ltime | log.Lshortfile)
	flag.Parse()

	if len(flag.Args()) != expectedNumArgs {
		flag.Usage()
		log.Fatalf("Error: expecting %d arguments, but got %d\n",
			expectedNumArgs, len(flag.Args()))
	}

	inFile := flag.Arg(0)
	outFile := flag.Arg(1)

	s := FilterSettings{
		InFile:     inFile,
		OutFile:    outFile,
		MatrixType: *matrixType,
		MinLength:  *minLength,
		MaxLength:  *maxLength,
	}

	pwmFilter(s)
}

func parseInfoArgs() {
	var expectedNumArgs int = 2
	var matrixType *string = flag.String("matrixType", "Weight", "Specify the type of position matrix. May be one of: Frequency, Probability, or Weight.")
	var pfmPseudocounts *float64 = flag.Float64("pfmPseudocounts", 0.1, "If a Position Frequency Matrix is provided, this pseudocount value will be applied when converting to a PWM.")
	flag.Usage = infoUsage
	log.SetFlags(log.Ldate | log.Ltime | log.Lshortfile)
	flag.Parse()

	if len(flag.Args()) != expectedNumArgs {
		flag.Usage()
		log.Fatalf("Error: expecting %d arguments, but got %d\n",
			expectedNumArgs, len(flag.Args()))
	}

	inFile := flag.Arg(0)
	outFile := flag.Arg(1)

	s := InfoSettings{
		InFile:       inFile,
		OutFile:      outFile,
		MatrixType:   *matrixType,
		Pseudocounts: *pfmPseudocounts,
	}

	pwmInfo(s)
}
