package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/interval"
	"log"
)

func bedEnrichments(method string, inFile string, secondFile string, noGapFile string, outFile string, verbose int, trimToRefGenome bool) {
	var err error

	if method != "exact" && method != "normalApproximate" && method != "upperBound" && method != "lowerBound" {
		log.Fatalf("Error: unknown method. Found: %s.", method)
	}

	if verbose > 0 {
		log.Println("Reading bed files.")
	}

	elementsOne := bed.ReadLite(inFile)
	elementsTwo := bed.ReadLite(secondFile)
	noGapRegions := bed.ReadLite(noGapFile)

	if trimToRefGenome {
		var overlap1, overlap2 []interval.Interval
		var trimmedE1 []bed.Bed = make([]bed.Bed, 0)
		var trimmedE2 []bed.Bed = make([]bed.Bed, 0)

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
			overlap1 = interval.Query(tree1, noGapRegions[i], "within")
			if overlap1 != nil {
				for j := range overlap1 {
					trimmedE1 = append(trimmedE1, overlap1[j].(bed.Bed))
				}
			}
			overlap2 = interval.Query(tree2, noGapRegions[i], "within")
			if overlap2 != nil {
				for j := range overlap2 {
					trimmedE2 = append(trimmedE2, overlap2[j].(bed.Bed))
				}
			}
		}
		elementsOne = trimmedE1
		elementsTwo = trimmedE2
	}

	if verbose > 0 {
		log.Println("Sorting bed files.")
	}

	bed.SortByCoord(elementsOne)
	bed.SortByCoord(elementsTwo)
	bed.SortByCoord(noGapRegions)

	if verbose > 0 {
		log.Println("Running preflight checks.")
	}

	//preflight checks: check for error in user input. Beds should not be self-overlapping.
	if bed.IsSelfOverlapping(elementsOne) {
		log.Fatalf("Elements in bedEnrichments must not be self-overlapping. Self-overlap found in inFile1.bed.")
	}
	if bed.IsSelfOverlapping(elementsTwo) {
		log.Fatalf("Elements in bedEnrichments must not be self-overlapping. Self-overlap found in inFile2.bed.")
	}
	if bed.IsSelfOverlapping(noGapRegions) {
		log.Fatalf("Elements in bedEnrichments must not be self-overlapping. Self-overlap found in noGap.bed.")
	}

	var summarySlice []float64
	overlapCount := bed.OverlapCount(elementsOne, elementsTwo)

	if verbose > 0 {
		log.Println("Calculating enrichment probabilities.")
	}

	switch method {
	case "exact":
		probs := bed.ElementOverlapProbabilities(elementsOne, elementsTwo, noGapRegions)
		summarySlice = bed.EnrichmentPValue(probs, overlapCount)
	case "normalApproximate":
		probs := bed.ElementOverlapProbabilities(elementsOne, elementsTwo, noGapRegions)
		summarySlice = bed.EnrichmentPValueApproximation(probs, overlapCount)
	case "upperBound":
		summarySlice = bed.EnrichmentPValueUpperBound(elementsOne, elementsTwo, noGapRegions, overlapCount, verbose)
	case "lowerBound":
		summarySlice = bed.EnrichmentPValueLowerBound(elementsOne, elementsTwo, noGapRegions, overlapCount, verbose)
	default:
		log.Fatalf("Error: unknown method. Found: %s.", method)
	}

	if verbose > 0 {
		log.Println("Done calculating enrichment. Writing to output.")
	}

	out := fileio.EasyCreate(outFile)
	_, err = fmt.Fprintf(out, "#Method\tFilename1\tFilename2\tLenElements1\tLenElements2\tOverlapCount\tDebugCheck\tExpectedOverlap\tEnrichment\tpValue\n")
	exception.PanicOnErr(err)
	_, err = fmt.Fprintf(out, "%s\t%s\t%s\t%d\t%d\t%d\t%f\t%f\t%f\t%e\n", method, inFile, secondFile, len(elementsOne), len(elementsTwo), overlapCount, summarySlice[0], summarySlice[1], float64(overlapCount)/summarySlice[1], summarySlice[2])
	exception.PanicOnErr(err)

	err = out.Close()
	exception.PanicOnErr(err)
}

func usage() {
	fmt.Print(
		"bedEnrichments - Returns the p-value of enrichment for overlaps between the elements in two input bed files.\n" +
			"noGap.bed represents a bed of all regions in the search space of the genome.\n" +
			"out.txt is in the form of a tab-separated value file with a header line starting with '#'.\n" +
			"Method specifies the type of calculation. Must match one of the following strings: 'exact', 'normalApproximate', 'upperBound', or 'lowerBound'\n" +
			"A brief explanation of each method of calculation is described starting on the next line.\n" +
			"exact: Calculates the exact p Value for enrichment between two bed files. NP-hard, computationally intractable for large datasets.\n" +
			"normalApproximate: (Recommended) Uses a normal approximation of the binomial distribution for the p Value. For large datasets, this method is rapid and should converge on the true p value, unless the p value is extremely small.\n" +
			"upperBound: Returns the most conservative exact estimate of the p value. Rapid, but the true p value is guaranteed to fall below this value.\n" +
			"lowerBound: Returns the lower bound exact estimate of the p value. Rapid, and the true p value is guaranteed to be above this value.\n" +
			"Usage:\n" +
			"bedEnrichments method elements1.bed elements2.bed noGap.bed out.txt\n")
	flag.PrintDefaults()
}

func main() {
	var expectedNumArgs int = 5

	var verbose *int = flag.Int("verbose", 0, "Set to 1 to reveal debug prints.")
	var trimToRefGenome *bool = flag.Bool("trimToRefGenome", false, "Ignores elements that do not lie within the reference genome, as defined by the noGap.bed file.")

	flag.Usage = usage
	//log.SetFlags(log.Ldate | log.Ltime | log.Lshortfile)
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

	bedEnrichments(method, inFile, secondFile, noGapFile, outFile, *verbose, *trimToRefGenome)
}
