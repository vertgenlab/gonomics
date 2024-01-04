package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/motif"
	"log"
	"os"
)

// FormatSettings defines the usage settings for the pwmTools format subcommand.
type FormatSettings struct {
	InFile      string
	OutFile     string
	InType      string
	OutType     string
	PseudoCount float64
}

// formatUsage defines the usage statement for the pwmTools format subcommand.
func formatUsage(formatFlags *flag.FlagSet) {
	fmt.Printf("pwmTools format - a tool for reformatting PFM/PPM/PWM files and for converting between these formats.\n" +
		"Usage:\n" +
		"pwmTools format in.pfm/ppm/pwm out.pfm/ppm/pwm\n" +
		"options:\n")
	formatFlags.PrintDefaults()
}

// parseFormatArgs is the main function of the pwmTools format subcommand. It parses options and runs the pwmFormat function.
func parseFormatArgs() {
	var expectedNumArgs int = 2
	var err error
	formatFlags := flag.NewFlagSet("format", flag.ExitOnError)
	var inType *string = formatFlags.String("inType", "Weight", "Specify the type of the input matrix file. Can be 'Frequency', 'Probability', or 'Weight'.")
	var outType *string = formatFlags.String("outType", "Frequency", "Specify the type of the output matrix file. Can be 'Frequency', 'Probability', or 'Weight'.")
	var pseudoCount *float64 = formatFlags.Float64("pfmPseudocounts", 0.1, "If a Position Frequency Matrix is provided, this pseudocount value will be applied when converting to a PWM or PPM.")
	err = formatFlags.Parse(os.Args[2:])
	exception.PanicOnErr(err)
	formatFlags.Usage = func() { formatUsage(formatFlags) }

	if len(formatFlags.Args()) != expectedNumArgs {
		formatFlags.Usage()
		log.Fatalf("Error: expecting %d arguments, but got %d\n",
			expectedNumArgs, len(formatFlags.Args()))
	}

	inFile := formatFlags.Arg(0)
	outFile := formatFlags.Arg(1)

	s := FormatSettings{
		InFile:      inFile,
		OutFile:     outFile,
		InType:      *inType,
		OutType:     *outType,
		PseudoCount: *pseudoCount,
	}

	pwmFormat(s)
}

// pwmFormat parses an input Position Matrix file and formats the file according to user-defined settings.
// Currently, this supports converting between PFM, PPM, and PWM motif representations.
func pwmFormat(s FormatSettings) {
	var records = motif.ReadJaspar(s.InFile, s.InType)
	switch s.InType {
	case "Frequency":
		switch s.OutType {
		case "Frequency":
			//nothing to do here, but also no need to fatal
		case "Probability":
			records = motif.PfmSliceToPpmSlice(records, s.PseudoCount)
		case "Weight":
			records = motif.PfmSliceToPpmSlice(records, s.PseudoCount)
			records = motif.PpmSliceToPwmSlice(records)
		default:
			log.Fatalf("Error: unrecognized output motif file format. Options are 'Frequency', 'Probability', and 'Weight'. Found: %v.\n", s.OutType)
		}
	case "Probability":
		switch s.OutType {
		case "Frequency":
			log.Fatalf("Error: Cannot convert a position probability matrix to a position frequency matrix.")
		case "Probability":
			//nothing to do here, but also no need to fatal
		case "Weight":
			records = motif.PpmSliceToPwmSlice(records)
		default:
			log.Fatalf("Error: Unrecognized output matrix type. Options are 'Frequency', 'Probability', or 'Weight'. Found: %v.\n", s.OutType)
		}
	case "Weight":
		switch s.OutFile {
		case "Frequency":
			log.Fatalf("Error: Cannot convert a position weight matrix to a position frequency matrix.")
		case "Probability":
			records = motif.PwmSliceToPpmSlice(records)
		case "Weight":
			//nothing to do here, but also no need to fatal
		default:
			log.Fatalf("Error: Unrecognized output matrix type. Options are 'Frequency', 'Probability', or 'Weight'. Found: %v.\n", s.OutType)
		}
	default:
		log.Fatalf("Error: unrecognized input matrix type. Options are 'Frequency', 'Probability', or 'Weight'. Found: %v.\n", s.InType)
	}

	motif.WriteJaspar(s.OutFile, records)
}
