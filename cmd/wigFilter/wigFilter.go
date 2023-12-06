// Command Group: "WIG Tools"

// Returns a filtered wig based on option parameters
package main

import (
	"flag"
	"fmt"
	"log"

	"github.com/vertgenlab/gonomics/wig"
)

type Settings struct {
	InFile       string
	ChromSizes   string
	OutFile      string
	Chrom        string
	DefaultValue float64
}

func wigFilter(s Settings) {
	var pass bool
	records := wig.ReadWholeGenome(s.InFile, s.ChromSizes, s.DefaultValue)
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
	wig.WriteMap(s.OutFile, answer)
}

func usage() {
	fmt.Print(
		"wigFilter - Returns a filtered wig based on option parameters.\n" +
			"Usage:\n" +
			"wigFilter input.wig chrom.sizes output.wig\n" +
			"options:\n")
	flag.PrintDefaults()
}

func main() {
	var expectedNumArgs int = 3
	var chrom *string = flag.String("chrom", "", "Retains wig entries with chromosome name matching this value.")
	var defaultValue *float64 = flag.Float64("defaultValue", 0, "Specifies the value in the wig for missing data.")

	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime | log.Lshortfile)
	flag.Parse()
	if len(flag.Args()) != expectedNumArgs {
		flag.Usage()
		log.Fatalf("Error: expecting %d arguments, but got %d\n",
			expectedNumArgs, len(flag.Args()))
	}

	inFile := flag.Arg(0)
	chromSizes := flag.Arg(1)
	outFile := flag.Arg(2)
	var s = Settings{
		InFile:       inFile,
		ChromSizes:   chromSizes,
		OutFile:      outFile,
		Chrom:        *chrom,
		DefaultValue: *defaultValue,
	}

	wigFilter(s)
}
