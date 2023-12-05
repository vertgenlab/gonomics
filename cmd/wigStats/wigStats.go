// Command Group: "WIG Tools"

// Provide coverage histogram for WIG format visualization files
package main

import (
	"flag"
	"fmt"
	"log"
	"math"

	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/wig"
)

type Settings struct {
	InFile           string
	ChromSizes       string
	NoGapFile        string
	OutFile          string
	MissingDataValue float64
}

func wigStats(s Settings) {
	var records []bed.Bed = bed.Read(s.NoGapFile)
	var w = wig.ReadWholeGenome(s.InFile, s.ChromSizes, s.MissingDataValue)
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

func usage() {
	fmt.Print(
		"wigStats - provide coverage histogram for WIG format visualization files.\n" +
			"wig values must be non-negative, unless the negative value is the missing data value.\n" +
			"Float values will be truncated to integers\n" +
			"Wigs must be fixed step of step size 1.\n" +
			"Usage:\n" +
			"wigStats in.wig chrom.sizes noGap.bed output.tsv\n" +
			"options:\n")
	flag.PrintDefaults()
}

func main() {
	var expectedNumArgs int = 4
	var missingDataValue *float64 = flag.Float64("missingDataValue", math.Inf(-1), "defines the value for missing data in the WIG file")
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
	noGapFile := flag.Arg(2)
	outFile := flag.Arg(3)

	s := Settings{
		InFile:           inFile,
		ChromSizes:       chromSizes,
		NoGapFile:        noGapFile,
		OutFile:          outFile,
		MissingDataValue: *missingDataValue,
	}

	wigStats(s)
}
