// Command Group: "Wig Tools"
//TODO make wigStats compatible with negative wig values

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
)

func wigStats(inFile string, noGapFile string, outFile string, missingDataValue float64) {
	var records []bed.Bed = bed.Read(noGapFile)
	var w []wig.Wig = wig.Read(inFile)
	var wigValues []float64
	var statValues, tmpValues []int
	var err error
	var i, j int
	statValues = make([]int, 1000)
	for i = range records {
		wigValues = wig.ChromToSlice(w, records[i].Chrom)
		for j = records[i].ChromStart; j < records[i].ChromEnd; j++ {
			if wigValues[j] == missingDataValue {
				continue
			}
			if int(wigValues[j]) > len(statValues) {
				tmpValues = make([]int, int(wigValues[j])-len(statValues)+1)
				statValues = append(statValues, tmpValues...) // concatenate two slices together
			}
			statValues[int(wigValues[j])]++
		}
	}
	out := fileio.EasyCreate(outFile)
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
		"wigStats - provide coverage statistics for WIG format visualization files.\n" +
			"Usage:\n" +
			"wigStats in.wig noGap.bed output.tsv\n" +
			"options:\n")
	flag.PrintDefaults()
}

func main() {
	var expectedNumArgs int = 3
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
	noGapFile := flag.Arg(1)
	outFile := flag.Arg(2)

	wigStats(inFile, noGapFile, outFile, *missingDataValue)
}
