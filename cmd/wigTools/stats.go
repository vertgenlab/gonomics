package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/wig"
	"log"
	"math"
	"os"
)

// StatsSettings defines the usage settings for the wigTools stats subcommand.
type StatsSettings struct {
	InFile           string
	ChromSizes       string
	NoGapFile        string
	OutFile          string
	MissingDataValue float64
}

// statsUsage defines the usage statement for the wigTools stats subcommand.
func statsUsage(statsFlags *flag.FlagSet) {
	fmt.Print(
		"wigTools stats - provide coverage histogram for WIG format visualization files.\n" +
			"wig values must be non-negative, unless the negative value is the missing data value.\n" +
			"Float values will be truncated to integers\n" +
			"Wigs must be fixed step of step size 1.\n" +
			"Usage:\n" +
			"wigTools stats in.wig chrom.sizes noGap.bed output.tsv\n" +
			"options:\n")
	statsFlags.PrintDefaults()
}

// parseStatsArgs is the main function of the wigTools stats subcommand. It parses options and runs the wigStats function.
func parseStatsArgs() {
	var expectedNumArgs int = 4
	var err error
	statsFlags := flag.NewFlagSet("stats", flag.ExitOnError)
	var missingDataValue *float64 = statsFlags.Float64("missingDataValue", math.Inf(-1), "defines the value for missing data in the WIG file")
	err = statsFlags.Parse(os.Args[2:])
	exception.PanicOnErr(err)
	statsFlags.Usage = func() { statsUsage(statsFlags) }

	if len(statsFlags.Args()) != expectedNumArgs {
		statsFlags.Usage()
		log.Fatalf("Error: expecting %d arguments, but got %d.\n", expectedNumArgs, len(statsFlags.Args()))
	}

	inFile := statsFlags.Arg(0)
	chromSizes := statsFlags.Arg(1)
	noGapFile := statsFlags.Arg(2)
	outFile := statsFlags.Arg(3)

	s := StatsSettings{
		InFile:           inFile,
		ChromSizes:       chromSizes,
		NoGapFile:        noGapFile,
		OutFile:          outFile,
		MissingDataValue: *missingDataValue,
	}
	wigStats(s)
}

func wigStats(s StatsSettings) {
	var records []bed.Bed = bed.Read(s.NoGapFile)
	var w = wig.Read(s.InFile, s.ChromSizes, s.MissingDataValue)
	var foundInMap bool
	var statValues, tmpValues []int
	var err error
	var i, j int
	statValues = make([]int, 1000)
	for i = range records {
		if _, foundInMap = w[records[i].Chrom]; !foundInMap {
			log.Fatalf("Error: chrom in bed entry: %s, not found in reference genome.\n", records[i].Chrom)
		}
		for j = records[i].ChromStart; j < records[i].ChromEnd; j++ {
			if w[records[i].Chrom].Values[j] == s.MissingDataValue {
				continue
			}
			if int(w[records[i].Chrom].Values[j]) >= len(statValues) {
				tmpValues = make([]int, int(w[records[i].Chrom].Values[j])-len(statValues)+1)
				statValues = append(statValues, tmpValues...) // concatenate two slices together
			}
			statValues[int(w[records[i].Chrom].Values[j])]++
		}
	}
	out := fileio.EasyCreate(s.OutFile)
	_, err = fmt.Fprintf(out, "coverage\tcount\n")
	exception.PanicOnErr(err)
	for i = range statValues {
		_, err = fmt.Fprintf(out, "%d\t%d\n", i, statValues[i])
		exception.PanicOnErr(err)
	}
	err = out.Close()
	exception.PanicOnErr(err)
}
