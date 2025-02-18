// Command Group: "Statistics & Population Genetics"

// Returns the p-value of enrichment and depletion for overlaps between the elements in two input files
package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/numbers"
	"log"

	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/interval"
	"github.com/vertgenlab/gonomics/interval/lift"
)

type Settings struct {
	Method            string
	InFile            string
	SecondFile        string
	SearchSpaceFile   string
	OutFile           string
	Verbose           int
	TrimToSearchSpace bool
	SecondFileList    string
	Relationship      string
}

func overlapEnrichments(s Settings) {
	var err error
	var searchSpaceIntervals []interval.Interval
	var searchSpaceTree map[string]*interval.IntervalNode
	if s.Method != "exact" && s.Method != "normalApproximate" && s.Method != "upperBound" && s.Method != "lowerBound" {
		log.Fatalf("Error: unknown method. Found: %s.", s.Method)
	}

	if s.Verbose > 0 {
		log.Println("Reading files.")
	}
	elementsOne := lift.GoRead(s.InFile)
	searchSpaceRegions := lift.GoRead(s.SearchSpaceFile)
	if s.TrimToSearchSpace {
		elementsOne = refGenomeTrim(elementsOne, searchSpaceRegions, s.Relationship)
	} else {
		// if TrimToSearchSpace is not active, throw an error if foreground elements fall outside background elements.
		searchSpaceIntervals = make([]interval.Interval, 0)
		for i := range searchSpaceRegions {
			searchSpaceIntervals = append(searchSpaceIntervals, searchSpaceRegions[i])
		}
		searchSpaceTree = interval.BuildTree(searchSpaceIntervals)
		for _, element := range elementsOne {
			if len(interval.Query(searchSpaceTree, element, "any")) == 0 {
				log.Fatalf("Error: foreground element from file 1 does not overlap search space. Please use 'trimToSearchSpace' to exclude this element. Offending element: %s.\n", element)
			}
		}
	}
	lift.SortByCoord(elementsOne)
	lift.SortByCoord(searchSpaceRegions)

	if s.Verbose > 0 {
		log.Println("Running preflight checks.")
	}

	if lift.IsSelfOverlapping(searchSpaceRegions, s.Verbose) {
		log.Fatalf("Elements in bedEnrichments must not be self-overlapping. Self-overlap found in %s.", s.SearchSpaceFile)
	}

	//preflight checks: check for error in user input. Beds should not be self-overlapping.
	if lift.IsSelfOverlapping(elementsOne, s.Verbose) {
		log.Fatalf("Elements in bedEnrichments must not be self-overlapping. Self-overlap found in %s.", s.InFile)
	}

	var secondFileList []string
	if s.SecondFileList == "" {
		secondFileList = []string{s.SecondFile}
	} else {
		secondFileList = fileio.Read(s.SecondFileList)
	}

	// Initialize output file
	out := fileio.EasyCreate(s.OutFile)
	_, err = fmt.Fprintf(out, "#Method\tFilename1\tFilename2\tLenElements1\tLenElements2\tOverlapCount\tDebugCheck\tExpectedOverlap\tEnrichment\tEnrichPValue\tDepletePValue\n")
	exception.PanicOnErr(err)

	for currSecondFile := range secondFileList {
		elementsTwo := lift.GoRead(secondFileList[currSecondFile])
		if s.TrimToSearchSpace {
			elementsTwo = refGenomeTrim(elementsTwo, searchSpaceRegions, s.Relationship)
		} else {
			for _, element := range elementsTwo {
				if len(interval.Query(searchSpaceTree, element, "any")) == 0 {
					log.Fatalf("Error: foreground element from file 2 does not overlap search space. Please use 'trimToSearchSpace' to exclude this element. Offending element: %s.\n", element)
				}
			}
		}
		lift.SortByCoord(elementsTwo)

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

		if s.Verbose > 0 {
			log.Println("Calculating enrichment probabilities.")
		}

		switch s.Method {
		case "exact":
			probs := lift.ElementOverlapProbabilities(elementsOne, elementsTwo, searchSpaceRegions)
			summarySlice = lift.EnrichmentPValueExact(probs, overlapCount)
		case "normalApproximate":
			probs := lift.ElementOverlapProbabilities(elementsOne, elementsTwo, searchSpaceRegions)
			summarySlice = lift.EnrichmentPValueApproximation(probs, overlapCount)
		case "upperBound":
			summarySlice = lift.EnrichmentPValueUpperBound(elementsOne, elementsTwo, searchSpaceRegions, overlapCount, s.Verbose)
		case "lowerBound":
			summarySlice = lift.EnrichmentPValueLowerBound(elementsOne, elementsTwo, searchSpaceRegions, overlapCount, s.Verbose)
		default:
			log.Fatalf("Error: unknown method. Found: %s.", s.Method)
		}

		if s.Verbose > 0 {
			log.Println("Done calculating enrichment. Writing to output.")
		}
		_, err = fmt.Fprintf(out, "%s\t%s\t%s\t%d\t%d\t%d\t%f\t%f\t%f\t%e\t%e\n", s.Method, s.InFile, secondFileList[currSecondFile], len(elementsOne), len(elementsTwo), overlapCount, summarySlice[0], summarySlice[1], float64(overlapCount)/summarySlice[1], summarySlice[2], summarySlice[3])
		exception.PanicOnErr(err)

	}
	err = out.Close()
	exception.PanicOnErr(err)
}

// refGenomeTrim trims elements in the input slice `unTrimmed` to only the portions that overlap with
// the regions specified in `noGapRegions`.
func refGenomeTrim(unTrimmed []lift.Lift, noGapRegions []lift.Lift, relationship string) []lift.Lift {
	var trimmed []lift.Lift = make([]lift.Lift, 0)
	var overlappingIntervals []interval.Interval
	var trimmedIntervals []interval.Interval

	var e1Intervals []interval.Interval
	for i := range unTrimmed {
		e1Intervals = append(e1Intervals, unTrimmed[i])
	}
	tree1 := interval.BuildTree(e1Intervals)

	for i := range noGapRegions {
		overlappingIntervals = interval.Query(tree1, noGapRegions[i], relationship)
		trimmedIntervals = trimIntervals(overlappingIntervals, noGapRegions[i])
		for _, t := range trimmedIntervals {
			trimmed = append(trimmed, t.(lift.Lift))
		}
	}
	return trimmed
}

// trimIntervals takes a set of intervals and trims each one to only the portion that overlaps with a given region.
func trimIntervals(overlappingIntervals []interval.Interval, region interval.Interval) []interval.Interval {
	var trimmed []interval.Interval
	for _, currInterval := range overlappingIntervals {
		trimmedStart := numbers.Max(currInterval.GetChromStart(), region.GetChromStart())
		trimmedEnd := numbers.Min(currInterval.GetChromEnd(), region.GetChromEnd())
		trimmedInterval := &bed.Bed{
			Chrom:             currInterval.GetChrom(),
			ChromStart:        trimmedStart,
			ChromEnd:          trimmedEnd,
			FieldsInitialized: 3,
		}
		trimmed = append(trimmed, trimmedInterval)
	}
	return trimmed
}

func usage() {
	fmt.Print(
		"overlapEnrichments - Returns the p-value of enrichment and depletion for overlaps between the elements in two input files.\n" +
			"search.lift represents a lift compatible file (current support for bed/vcf) of all regions in the search space of the genome. This is often a noGap.bed for whole genome backgrounds.\n" +
			"Genomic elements are expected to lie within the search space, and an error will occur if elements outside the search space are detected. Alternatively, the user can specify 'trimToSearchSpace'\n" +
			"to ignore elements outside of the search space.\n" +
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
			"overlapEnrichments method elements1.lift elements2.lift searchSpace.lift out.txt\n")
	flag.PrintDefaults()
}

func main() {
	var expectedNumArgs int = 5
	var verbose *int = flag.Int("verbose", 0, "Set to 1 to reveal debug prints.")
	var trimToSearchSpace *bool = flag.Bool("trimToSearchSpace", false, "Ignores elements that do not lie within the search space, as defined by the searchSpace.lift file.")
	var secondFileList *string = flag.String("secondFileList", "", "Specify a list of query files to calculate enrichments against the first file. Note that while using this option the command will ignore the elements2.lift argument.")
	var relationship *string = flag.String("relationship", "within", "Specify an overlap relationship for the trimToSearchSpace option. 'all' is more permissive than the default 'within'.")

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
	searchSpaceFile := flag.Arg(3)
	outFile := flag.Arg(4)

	s := Settings{
		Method:            method,
		InFile:            inFile,
		SecondFile:        secondFile,
		SearchSpaceFile:   searchSpaceFile,
		OutFile:           outFile,
		Verbose:           *verbose,
		TrimToSearchSpace: *trimToSearchSpace,
		SecondFileList:    *secondFileList,
		Relationship:      *relationship,
	}

	overlapEnrichments(s)
}
