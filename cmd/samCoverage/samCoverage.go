// Command Group: "SAM Tools"

// Calculates genome coverage as the quotient of aligned bases in a sequencing dataset to the total length of ungapped genomic regions in the reference genome
package main

import (
	"flag"
	"fmt"
	"log"

	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/numbers"
	"github.com/vertgenlab/gonomics/numbers/fit"
	"github.com/vertgenlab/gonomics/sam"
)

func samCoverage(samFileName string, outFile string, countNinDepth bool) {
	out := fileio.EasyCreate(outFile)
	fmt.Fprintf(out, "Coverage\tPileups\tGroup\n")
	data, header := sam.GoReadToChan(samFileName)
	piles := sam.GoPileup(data, header, false, nil, nil)
	histogram := make([]float64, 30)
	observations := make([]float64, 1000)
	totalCount := 0

	for p := range piles {
		//if else structure to calculate pile depth
		depth := TotalDepth(p, countNinDepth)
		if depth >= len(histogram) {
			// extend histogram to the length of depth
			histogramBuffer := make([]float64, depth+10) // this is a slice of 0s used to extend histogram when required
			copy(histogramBuffer, histogram)
			histogram = histogramBuffer
		}
		if totalCount >= len(observations) {
			observationsBuffer := make([]float64, len(observations)+1000)
			copy(observationsBuffer, observations)
			observations = observationsBuffer
		}
		histogram[depth]++
		observations[totalCount] = float64(depth)
		totalCount++
	}

	lambda := fit.Poisson(observations)
	for i, pileups := range histogram {
		fmt.Fprintf(out, "%v\t%v\tEmpirical\n", i, pileups)
		y := numbers.PoissonDist(i, lambda)
		fmt.Println(y)
		fmt.Fprintf(out, "%v\t%v\tExpected\n", i, y*float64(totalCount))
	}
	fmt.Printf("lambda:%v\n", lambda)
	out.Close()
}

func TotalDepth(p sam.Pile, countNinDepth bool) int {
	depth := p.CountF[dna.A] + p.CountF[dna.C] + p.CountF[dna.G] + p.CountF[dna.T] + p.CountR[dna.A] + p.CountR[dna.C] + p.CountR[dna.G] + p.CountR[dna.T]
	if countNinDepth {
		depth += p.CountF[dna.N] + p.CountR[dna.N]
	}
	return depth
}

func usage() {
	fmt.Print(
		"samCoverage - Calculates genome coverage as the quotient of aligned bases in a sequencing dataset to the total length of ungapped genomic regions in the reference genome.\n" +
			"Usage:\n" +
			"samCoverage input.sam nogap.bed outfile.txt\n" +
			"options:\n")
	flag.PrintDefaults()
}

func main() {
	var countNinDepth *bool = flag.Bool("countNinDepth", true, "If true, count 'N' reads towards total depth of pileups.")
	var expectedNumArgs int = 3
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

	samCoverage(samFileName, outFile, *countNinDepth)
}
