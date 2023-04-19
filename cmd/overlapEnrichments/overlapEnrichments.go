package main

import (
	"flag"
	"fmt"
	"log"

	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/interval"
	"github.com/vertgenlab/gonomics/interval/lift"
)

func overlapEnrichments(method string, inFile string, secondFile string, noGapFile string, outFile string, verbose int, trimToRefGenome bool) {
	var err error

	if method != "exact" && method != "normalApproximate" && method != "upperBound" && method != "lowerBound" {
		log.Fatalf("Error: unknown method. Found: %s.", method)
	}

	if verbose > 0 {
		log.Println("Reading files.")
	}

	elementsOne := lift.GoRead(inFile)
	elementsTwo := lift.GoRead(secondFile)
	noGapRegions := lift.GoRead(noGapFile)

	if trimToRefGenome {
		var overlap1, overlap2 []lift.Lift
		var trimmedE1 []lift.Lift = make([]lift.Lift, 0)
		var trimmedE2 []lift.Lift = make([]lift.Lift, 0)

		var e1Intervals []interval.Interval
		for i := range elementsOne {
			e1Intervals = append(e1Intervals, elementsOne[i])
		}
		tree1 := interval.BuildTree(e1Intervals)
		var e2Intervals []interval.Interval
		for i := range elementsTwo {
			e2Intervals = append(e2Intervals, elementsTwo[i])
		}
		tree2 := interval.BuildTree(e2Intervals)

		for i := range noGapRegions {
			overlap1 = lift.IntervalSliceToLift(interval.Query(tree1, noGapRegions[i], "within"))
			for j := range overlap1 {
				trimmedE1 = append(trimmedE1, overlap1[j])
			}
			overlap2 = lift.IntervalSliceToLift(interval.Query(tree2, noGapRegions[i], "within"))
			for j := range overlap2 {
				trimmedE2 = append(trimmedE2, overlap2[j])
			}
		}
		elementsOne = trimmedE1
		elementsTwo = trimmedE2
	}

	if verbose > 0 {
		log.Println("Sorting bed files.")
	}

	lift.SortByCoord(elementsOne)
	lift.SortByCoord(elementsTwo)
	lift.SortByCoord(noGapRegions)

	if verbose > 0 {
		log.Println("Running preflight checks.")
	}

	if lift.IsSelfOverlapping(noGapRegions, verbose) {
		log.Fatalf("Elements in bedEnrichments must not be self-overlapping. Self-overlap found in %s.", noGapFile)
	}

	//preflight checks: check for error in user input. Beds should not be self-overlapping.
	if lift.IsSelfOverlapping(elementsOne, verbose) {
		log.Fatalf("Elements in bedEnrichments must not be self-overlapping. Self-overlap found in %s.", inFile)
	}

	// Check below removed to enable overlapping elements for second bed region.
	// The first file is viewed as elements/annotations of the human genome that are labels on the genomic coordinates.
	// The idea of a region being labeled multiple times does not fit cleanly into the math because it is a little vague
	// how a single element from list two overlapping multiple elements from list 1 should be counted.
	// Right now we count it as a single "success" where the maximum number of "successes" is the number of elements in file2.
	// File2 really represents sampling from this genome that has been annotated with file1, so there it is easier
	// to have them overlapping, since they are treated as independent draws
	//
	//if lift.IsSelfOverlapping(elementsTwo, verbose) {
	//	log.Fatalf("Elements in bedEnrichments must not be self-overlapping. Self-overlap found in %s.", secondFile)
	//}

	var summarySlice []float64
	overlapCount := lift.OverlapCount(elementsTwo, elementsOne)

	if verbose > 0 {
		log.Println("Calculating enrichment probabilities.")
	}

	switch method {
	case "exact":
		probs := lift.ElementOverlapProbabilities(elementsOne, elementsTwo, noGapRegions)
		summarySlice = lift.EnrichmentPValueExact(probs, overlapCount)
	case "normalApproximate":
		probs := lift.ElementOverlapProbabilities(elementsOne, elementsTwo, noGapRegions)
		summarySlice = lift.EnrichmentPValueApproximation(probs, overlapCount)
	case "upperBound":
		summarySlice = lift.EnrichmentPValueUpperBound(elementsOne, elementsTwo, noGapRegions, overlapCount, verbose)
	case "lowerBound":
		summarySlice = lift.EnrichmentPValueLowerBound(elementsOne, elementsTwo, noGapRegions, overlapCount, verbose)
	default:
		log.Fatalf("Error: unknown method. Found: %s.", method)
	}

	if verbose > 0 {
		log.Println("Done calculating enrichment. Writing to output.")
	}

	out := fileio.EasyCreate(outFile)
	_, err = fmt.Fprintf(out, "#Method\tFilename1\tFilename2\tLenElements1\tLenElements2\tOverlapCount\tDebugCheck\tExpectedOverlap\tEnrichment\tEnrichPValue\tDepletePValue\n")
	exception.PanicOnErr(err)
	_, err = fmt.Fprintf(out, "%s\t%s\t%s\t%d\t%d\t%d\t%f\t%f\t%f\t%e\t%e\n", method, inFile, secondFile, len(elementsOne), len(elementsTwo), overlapCount, summarySlice[0], summarySlice[1], float64(overlapCount)/summarySlice[1], summarySlice[2], summarySlice[3])
	exception.PanicOnErr(err)

	err = out.Close()
	exception.PanicOnErr(err)
}

func usage() {
	fmt.Print(
		"overlapEnrichments - Returns the p-value of enrichment and depletion for overlaps between the elements in two input files.\n" +
			"noGap.lift represents a lift compatible file (current support for bed/vcf) of all regions in the search space of the genome.\n" +
			"out.txt is in the form of a tab-separated value file with a header line starting with '#'.\n" +
			"Calculates enrichment of the number of elements in set 2 that have any overlap in set 1.\n" +
			"Number of overlaps reported is the number of elements in set 2 that have any overlap with set 1. This will be asymmetric if sets one and two are swapped as arguments.\n" +
			"Method specifies the type of calculation. Must match one of the following strings: 'exact', 'normalApproximate', 'upperBound', or 'lowerBound'\n" +
			"A brief explanation of each method of calculation is described starting on the next line.\n" +
			"exact: Calculates the exact p-Value for enrichment/depletion between two files. NP-hard, computationally intractable for large datasets.\n" +
			"normalApproximate: (Recommended) Uses a normal approximation of the binomial distribution for the p-Value. For large datasets, this method is rapid and should converge on the true p value, unless the p value is extremely small.\n" +
			"upperBound: Returns the most conservative exact estimate of the p-Value. Rapid, but the true p-Value is guaranteed to fall below this value.\n" +
			"lowerBound: Returns the lower bound exact estimate of the p-Value. Rapid, and the true p value is guaranteed to be above this value.\n" +
			"Usage:\n" +
			"overlapEnrichments method elements1.lift elements2.lift noGap.lift out.txt\n")
	flag.PrintDefaults()
}

func main() {
	var expectedNumArgs int = 5

	var verbose *int = flag.Int("verbose", 0, "Set to 1 to reveal debug prints.")
	var trimToRefGenome *bool = flag.Bool("trimToRefGenome", false, "Ignores elements that do not lie within the reference genome, as defined by the noGap.bed file.")

	flag.Usage = usage
	log.SetFlags(0)
	flag.Parse()

	if len(flag.Args()) != expectedNumArgs {
		flag.Usage()
		log.Fatalf("Error: expecting %d arguments, but got %d\n",
			expectedNumArgs, len(flag.Args()))
	}
	method := flag.Arg(0)
	inFile := flag.Arg(1)
	secondFile := flag.Arg(2)
	noGapFile := flag.Arg(3)
	outFile := flag.Arg(4)

	overlapEnrichments(method, inFile, secondFile, noGapFile, outFile, *verbose, *trimToRefGenome)
}
