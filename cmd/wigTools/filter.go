package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/wig"
	"log"
	"os"
)

// FilterSettings defines the usage settings for the wigTools filter subcommand.
type FilterSettings struct {
	InFile       string
	ChromSizes   string
	OutFile      string
	Chrom        string
	DefaultValue float64
}

func filterUsage(filterFlags *flag.FlagSet) {
	fmt.Printf("wigTools filter - returns a filtered wig based on user-specified parameters.\n" +
		"Usage:\n" +
		"wigTools filter input.wig chrom.sizes output.wig\n" +
		"options:\n")
	filterFlags.PrintDefaults()
}

func parseFilterArgs() {
	var expectedNumArgs int = 3
	var err error
	filterFlags := flag.NewFlagSet("filter", flag.ExitOnError)
	var chrom *string = filterFlags.String("chrom", "", "Retains wig entries with chromosome name matching this value.")
	var defaultValue *float64 = filterFlags.Float64("defaultValue", 0, "Specifies the value in the wig for missing data.")
	err = filterFlags.Parse(os.Args[2:])
	exception.PanicOnErr(err)
	filterFlags.Usage = func() { filterUsage(filterFlags) }
	if len(filterFlags.Args()) != expectedNumArgs {
		filterFlags.Usage()
		log.Fatalf("Error: expecting %d arguments, but got %d.\n", expectedNumArgs, len(filterFlags.Args()))
	}

	inFile := filterFlags.Arg(0)
	sizes := filterFlags.Arg(1)
	outFile := filterFlags.Arg(2)

	s := FilterSettings{
		InFile:       inFile,
		ChromSizes:   sizes,
		OutFile:      outFile,
		Chrom:        *chrom,
		DefaultValue: *defaultValue,
	}

	wigFilter(s)
}

func wigFilter(s FilterSettings) {
	var pass bool
	records := wig.Read(s.InFile, s.ChromSizes, s.DefaultValue)
	var answer = make(map[string]wig.Wig)

	for currKey := range records {
		pass = true
		if s.Chrom != "" && records[currKey].Chrom != s.Chrom {
			pass = false
		}
		if pass {
			answer[currKey] = records[currKey]
		}
	}
	wig.Write(s.OutFile, answer)
}
