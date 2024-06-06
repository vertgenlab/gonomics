// Command Group: "SAM Tools"

// Generates a count histogram of genome coverage.
package main

import (
	"flag"
	"fmt"
	"log"

	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/numbers"
	"github.com/vertgenlab/gonomics/numbers/fit"
	"github.com/vertgenlab/gonomics/sam"
)

type Settings struct {
	SamFileName   string
	OutFile       string
	CountNinDepth bool
	Verbose       int
}

// samCoverage generates a count histogram of genome coverage.
func samCoverage(s Settings) {
	out := fileio.EasyCreate(s.OutFile)
	_, err := fmt.Fprintf(out, "Coverage\tPileups\tGroup\n")
	exception.PanicOnErr(err)
	data, header := sam.GoReadToChan(s.SamFileName)
	piles := sam.GoPileup(data, header, false, nil, nil)
	histogram := make([]int, 30)
	totalCount := 0

	for p := range piles {
		depth := TotalDepth(p, s.CountNinDepth)
		if depth >= len(histogram) {
			// extend histogram to the length of depth
			histogramBuffer := make([]int, depth+10) // this is a slice of 0s used to extend histogram when required
			copy(histogramBuffer, histogram)
			histogram = histogramBuffer
		}
		histogram[depth]++
		totalCount++
	}
	lambda := fit.PoissonHistogram(histogram)
	for i, pileups := range histogram {
		_, err = fmt.Fprintf(out, "%v\t%v\tEmpirical\n", i, pileups)
		exception.PanicOnErr(err)
		y, outlier := numbers.PoissonDist(i, lambda, false)
		if !outlier {
			_, err = fmt.Fprintf(out, "%v\t%.6g\tExpected\n", i, y*float64(totalCount))
			exception.PanicOnErr(err)
		}
	}
	if s.Verbose > 0 {
		log.Printf("lambda:%v\n", lambda)
	}
	out.Close()
}

// TotalDepth adds up the counts for all the bases at the given pileup and returns the total count as the depth
func TotalDepth(p sam.Pile, countNinDepth bool) int {
	depth := p.CountF[dna.A] + p.CountF[dna.C] + p.CountF[dna.G] + p.CountF[dna.T] + p.CountR[dna.A] + p.CountR[dna.C] + p.CountR[dna.G] + p.CountR[dna.T]
	if countNinDepth {
		depth += p.CountF[dna.N] + p.CountR[dna.N]
	}
	return depth
}

func usage() {
	fmt.Print(
		"samCoverage - Generates a count histogram of genome coverage.\n" +
			"Usage:\n" +
			"samCoverage input.sam outfile.txt\n" +
			"options:\n")
	flag.PrintDefaults()
}

func main() {
	var countNinDepth *bool = flag.Bool("countNinDepth", true, "If true, count 'N' reads towards total depth of pileups.")
	var verbose *int = flag.Int("verbose", 0, "Set to 1 to reveal debug prints. Verbose in this program reports Poisson parameter lambda.")
	var expectedNumArgs int = 2
	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime | log.Lshortfile)
	flag.Parse()

	if len(flag.Args()) != expectedNumArgs {
		flag.Usage()
		log.Fatalf("Error: expecting %d arguments, but got %d\n",
			expectedNumArgs, len(flag.Args()))
	}

	samFileName := flag.Arg(0)
	outFile := flag.Arg(1)

	s := Settings{
		SamFileName:   samFileName,
		OutFile:       outFile,
		CountNinDepth: *countNinDepth,
		Verbose:       *verbose,
	}

	samCoverage(s)
}
