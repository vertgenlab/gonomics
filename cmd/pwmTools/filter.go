package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/motif"
	"github.com/vertgenlab/gonomics/numbers"
	"log"
	"os"
)

// FilterSettings defines the usage settings for the pwmTools filter subcommand.
type FilterSettings struct {
	InFile     string
	OutFile    string
	MatrixType string
	MinLength  int
	MaxLength  int
}

// filterUsage defines the usage statement for the pwmTools filter subcommand.
func filterUsage(filterFlags *flag.FlagSet) {
	fmt.Printf("pwmTools filter - a tool for filtering PWM records.\n" +
		"Usage:\n" +
		"pwmTools filter in.pwm out.pwm\n" +
		"options:\n")
	filterFlags.PrintDefaults()
}

// parseFilterArgs is the main function for the pwmTools filter subcommand. It parses options and to run the function pwmFilter.
func parseFilterArgs() {
	var expectedNumArgs int = 2
	var err error
	filterFlags := flag.NewFlagSet("filter", flag.ExitOnError)
	var minLength *int = filterFlags.Int("minLength", 0, "Specify the minimum length of PWMs to be retained in the output.")
	var maxLength *int = filterFlags.Int("maxLength", numbers.MaxInt, "Specify the maximum length of PWMs to be retained in the output.")
	var matrixType *string = filterFlags.String("matrixType", "Frequency", "Specify the type of position matrix. May be one of the following: Frequency, Probability, Weight.")
	err = filterFlags.Parse(os.Args[2:])
	exception.PanicOnErr(err)
	filterFlags.Usage = func() { filterUsage(filterFlags) }

	if len(filterFlags.Args()) != expectedNumArgs {
		filterFlags.Usage()
		log.Fatalf("Error: expecting %d arguments, but got %d\n",
			expectedNumArgs, len(filterFlags.Args()))
	}

	inFile := filterFlags.Arg(0)
	outFile := filterFlags.Arg(1)

	s := FilterSettings{
		InFile:     inFile,
		OutFile:    outFile,
		MatrixType: *matrixType,
		MinLength:  *minLength,
		MaxLength:  *maxLength,
	}

	pwmFilter(s)
}

// pwmFilter parses an input PositionMatrix file and retains entries in an output file that pass filter criteria.
func pwmFilter(s FilterSettings) {
	var err error
	var pass bool
	records := motif.ReadJaspar(s.InFile, s.MatrixType)
	out := fileio.EasyCreate(s.OutFile)

	for i := range records {
		pass = true
		if len(records[i].Mat[0]) < s.MinLength {
			pass = false
		} else if len(records[i].Mat[0]) > s.MaxLength {
			pass = false
		}
		if pass {
			motif.WritePositionMatrixJaspar(out, records[i])
		}
	}

	err = out.Close()
	exception.PanicOnErr(err)
}
