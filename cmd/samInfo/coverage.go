package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/numbers"
	"github.com/vertgenlab/gonomics/numbers/fit"
	"github.com/vertgenlab/gonomics/sam"
	"log"
	"os"
)

type CoverageSettings struct {
	SamFileName      string
	HistogramOutFile string
	StatSummaryFile  string
	HighEndFilter    float64
	CountNinDepth    bool
	Verbose          int
}

func coverageUsage(coverageFlags *flag.FlagSet) {
	fmt.Print(
		"samInfo coverage - Generates a count histogram of genome coverage.\n" +
			"Usage:\n" +
			"samInfo coverage input.sam histogram.txt statSummary.txt\n" +
			"options:\n")
	coverageFlags.PrintDefaults()
}

func parseCoverageArgs() {
	var expectedNumArgs int = 3
	var err error
	coverageFlags := flag.NewFlagSet("coverage", flag.ExitOnError)
	err = coverageFlags.Parse(os.Args[2:])
	exception.PanicOnErr(err)
	coverageFlags.Usage = func() { coverageUsage(coverageFlags) }
	if len(coverageFlags.Args()) != expectedNumArgs {
		coverageFlags.Usage()
		log.Fatalf("Error: expecting %d arguments, but got %d\n",
			expectedNumArgs, len(coverageFlags.Args()))
	}

	var countNinDepth *bool = coverageFlags.Bool("countNinDepth", true, "If true, count 'N' reads towards total depth of pileups.")
	var verbose *int = coverageFlags.Int("verbose", 0, "Set to 1 to reveal debug prints. Verbose in this program reports Poisson parameter lambda.")
	var highEndFilter *float64 = coverageFlags.Float64("highEndFilter", 0.001, "Percent threshold from right end of distribution to be filtered out")
	samFileName := flag.Arg(0)
	outHistFile := flag.Arg(1)
	outStatFile := flag.Arg(2)
	s := CoverageSettings{
		SamFileName:      samFileName,
		HistogramOutFile: outHistFile,
		StatSummaryFile:  outStatFile,
		HighEndFilter:    *highEndFilter,
		CountNinDepth:    *countNinDepth,
		Verbose:          *verbose,
	}
	samCoverage(s)
}

// TotalDepth adds up the counts for all the bases at the given pileup and returns the total count as the depth
func TotalDepth(p sam.Pile, countNinDepth bool) int {
	depth := p.CountF[dna.A] + p.CountF[dna.C] + p.CountF[dna.G] + p.CountF[dna.T] + p.CountR[dna.A] + p.CountR[dna.C] + p.CountR[dna.G] + p.CountR[dna.T]
	if countNinDepth {
		depth += p.CountF[dna.N] + p.CountR[dna.N]
	}
	return depth
}

// ThresholdCalc determines the coverage value associated with the filter threshold as a proportion (e.g. 0.001)
func ThresholdCalc(threshold float64, covHistogram []int, totalCount float64) int {
	targetNum := totalCount - totalCount*threshold
	observations := totalCount
	index := int(len(covHistogram) - 1)
	for observations > targetNum {
		observations -= float64(covHistogram[index])
		index--
	}
	return index
}

func samCoverage(s CoverageSettings) {
	outHist := fileio.EasyCreate(s.HistogramOutFile)
	statSummary := fileio.EasyCreate(s.StatSummaryFile)
	_, err := fmt.Fprintf(outHist, "Coverage\tPileups\tGroup\tFilename\n")
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
	coverageThreshold := ThresholdCalc(s.HighEndFilter, histogram, float64(totalCount))
	_, err = fmt.Fprintf(statSummary, "Lambda\t%v\nCoverageThreshold\t%v\n", lambda, coverageThreshold)
	exception.PanicOnErr(err)
	for i, pileups := range histogram {
		_, err = fmt.Fprintf(outHist, "%v\t%v\tEmpirical\t%v\n", i, pileups, s.SamFileName)
		exception.PanicOnErr(err)
		y, outlier := numbers.PoissonDist(i, lambda, false)
		if !outlier {
			_, err = fmt.Fprintf(outHist, "%v\t%.6g\tExpected\t%v\n", i, y*float64(totalCount), s.SamFileName)
			exception.PanicOnErr(err)
		}
	}
	if s.Verbose > 0 {
		log.Printf("Lambda: %v\n", lambda)
		log.Printf("Coverage Threshold for High-End Filter: %v\n", coverageThreshold)
	}
	err = outHist.Close()
	exception.PanicOnErr(err)
	err = statSummary.Close()
	exception.PanicOnErr(err)
}
