// Command Group: "Data Conversion"

// Convert HiC contact maps in straw format to bedpe contact peak calls
package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/bed/bedpe"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/hic"
	"github.com/vertgenlab/gonomics/numbers"
	"github.com/vertgenlab/gonomics/numbers/fit"
	"log"
	"math"
	"strings"
)

type Settings struct {
	FileList         string
	OutFile          string
	BinSize          int
	RStart           float64
	PStart           float64
	RStep            float64
	PStep            float64
	MinCutoff        int
	MinBinDistance   int
	Fdr              float64
	ContactScoreFile string
	FitStatsFile     string
}

func strawToBedpe(s Settings) {
	var currStrawChan <-chan hic.Straw
	var words []string
	var binDistance int
	var currDistance float64
	var contactScoreCache [][]int = make([][]int, 2)
	contactScoreCache[0] = make([]int, 1)
	contactScoreCache[1] = make([]int, 1)
	var tmpContactScoreCache [][]int
	var tmpColumn []int
	var err error
	var firstTime, foundInMap bool
	var currChrom string
	var searchSpaceMins = make(map[string]int)
	var searchSpaceMaxes = make(map[string]int)
	var currBedPe bedpe.BedPe

	lines := fileio.Read(s.FileList)
	for i := range lines {
		words = strings.Split(lines[i], "\t")
		currStrawChan = hic.GoReadToChan(words[0])
		currChrom = words[1]
		firstTime = true
		for currStraw := range currStrawChan {
			if firstTime {
				if _, foundInMap = searchSpaceMins[currChrom]; foundInMap {
					searchSpaceMins[currChrom] = numbers.Min(searchSpaceMins[currChrom], numbers.Min(currStraw.Bin1Start, currStraw.Bin2Start))
					searchSpaceMaxes[currChrom] = numbers.Max(searchSpaceMaxes[currChrom], numbers.Max(currStraw.Bin1Start, currStraw.Bin2Start))
				} else {
					searchSpaceMins[currChrom] = numbers.Min(currStraw.Bin1Start, currStraw.Bin2Start)
					searchSpaceMaxes[currChrom] = numbers.Max(currStraw.Bin1Start, currStraw.Bin2Start)
				}
				firstTime = false
			} else {
				if currStraw.Bin1Start < searchSpaceMins[currChrom] {
					searchSpaceMins[currChrom] = currStraw.Bin1Start
				}
				if currStraw.Bin2Start < searchSpaceMins[currChrom] {
					searchSpaceMins[currChrom] = currStraw.Bin2Start
				}
				if currStraw.Bin1Start > searchSpaceMaxes[currChrom] {
					searchSpaceMaxes[currChrom] = currStraw.Bin1Start
				}
				if currStraw.Bin2Start > searchSpaceMaxes[currChrom] {
					searchSpaceMaxes[currChrom] = currStraw.Bin2Start
				}
			}

			currDistance = math.Abs(float64(currStraw.Bin1Start) - float64(currStraw.Bin2Start))
			if int(currDistance)%s.BinSize != 0 {
				log.Fatalf("Error: Distance between two straw ends: %v is not a multiple of the bin size: %v.\n", currDistance, s.BinSize)
			}
			binDistance = int(currDistance) / s.BinSize
			// if we need more room in the cache rows, extend the cache
			if binDistance > len(contactScoreCache)-1 {
				tmpContactScoreCache = make([][]int, binDistance+1)
				copy(tmpContactScoreCache, contactScoreCache)
				contactScoreCache = tmpContactScoreCache
			}

			if contactScoreCache[binDistance] == nil {
				contactScoreCache[binDistance] = make([]int, 1)
			}
			// if we need more room in the cache columns, extend the cache
			if currStraw.ContactScore > len(contactScoreCache[binDistance])-1 {
				tmpColumn = make([]int, currStraw.ContactScore+1)
				copy(tmpColumn, contactScoreCache[binDistance])
				contactScoreCache[binDistance] = tmpColumn
			}
			contactScoreCache[binDistance][currStraw.ContactScore]++
		}
	}

	if s.ContactScoreFile != "" {
		printContactScoreCacheToFile(contactScoreCache, s)
	}

	// this cache stores in comparisonCountCache[i] the total number of bins at each bin distance i.
	// Used for FDR correction as the number of comparisons.
	comparisonCountCache := makeComparisonCountCache(contactScoreCache, searchSpaceMins, searchSpaceMaxes, s)
	cutoffCache := calculateBenjamaniHochbergCutoff(contactScoreCache, s, comparisonCountCache)

	out := fileio.EasyCreate(s.OutFile)
	// now we read all the straw files again and filter based on our cutoff.
	for i := range lines {
		words = strings.Split(lines[i], "\t")
		currStrawChan = hic.GoReadToChan(words[0])
		currChrom = words[1]
		for currStraw := range currStrawChan {
			currDistance = math.Abs(float64(currStraw.Bin1Start) - float64(currStraw.Bin2Start))
			if int(currDistance)%s.BinSize != 0 {
				log.Fatalf("Error: Distance between two straw ends: %v is not a multiple of the bin size: %v.\n", currDistance, s.BinSize)
			}
			binDistance = int(currDistance) / s.BinSize
			if binDistance >= s.MinBinDistance && currStraw.ContactScore > cutoffCache[binDistance] {
				currBedPe = bedpe.BedPe{
					A: bed.Bed{Chrom: currChrom,
						ChromStart:        currStraw.Bin1Start,
						ChromEnd:          currStraw.Bin1Start + s.BinSize,
						Score:             currStraw.ContactScore,
						FieldsInitialized: 8,
					},
					B: bed.Bed{Chrom: currChrom,
						ChromStart:        currStraw.Bin2Start,
						ChromEnd:          currStraw.Bin2Start + s.BinSize,
						Score:             currStraw.ContactScore,
						FieldsInitialized: 8,
					},
				}
				bedpe.WriteToFileHandle(out, currBedPe)
			}
		}
	}

	err = out.Close()
	exception.PanicOnErr(err)
}

func printContactScoreCacheToFile(contactScoreCache [][]int, s Settings) {
	var err error
	out := fileio.EasyCreate(s.ContactScoreFile)
	var currScore, currBinDistance int
	headerString := "Score"
	for i := range contactScoreCache {
		headerString += fmt.Sprintf("\tBinDistance.%v.Count", i)
	}
	_, err = fmt.Fprintf(out, "%v\n", headerString)
	exception.PanicOnErr(err)
	for currScore = range contactScoreCache[0] {
		_, err = fmt.Fprintf(out, "%v", currScore)
		for currBinDistance = range contactScoreCache {
			if currScore < len(contactScoreCache[currBinDistance]) {
				_, err = fmt.Fprintf(out, "\t%v", contactScoreCache[currBinDistance][currScore])
				exception.PanicOnErr(err)
			} else {
				_, err = fmt.Fprintf(out, "\t0")
				exception.PanicOnErr(err)
			}

		}
		_, err = fmt.Fprintf(out, "\n")
		exception.PanicOnErr(err)
	}

	err = out.Close()
	exception.PanicOnErr(err)
}

// makecomparisonCountCache is a helper function that calculates the number of comparisons made for each binDistance.
func makeComparisonCountCache(contactScoreCache [][]int, searchSpaceMins map[string]int, searchSpaceMaxes map[string]int, s Settings) []int {
	var totalWindows int
	var currKey string
	var comparisonCountCache []int = make([]int, len(contactScoreCache))
	for i := range comparisonCountCache {
		totalWindows = 0
		for currKey, _ = range searchSpaceMins {
			//fmt.Printf("currKey: %v. searchSpaceMins[currKey]: %v. searchSpaceMaxes[currKey]: %v. BinDistance: %v.\n", currKey, searchSpaceMins[currKey], searchSpaceMaxes[currKey], i)
			// this calculation is the number of valid contact windows on a chromosome. ChromSize/binSize is windows
			// on chromosome, minus i, the binDistance.
			totalWindows += (searchSpaceMaxes[currKey]-searchSpaceMins[currKey])/s.BinSize - i
		}
		comparisonCountCache[i] = totalWindows
	}
	return comparisonCountCache
}

// calculateBenjamaniHochbergCutoff takes in the contactScoreCache, and number of comparisons defined in the comparisonCountCache
// to return the cutoff point for each binDistance where the FDR-adjusted p value falls below the user-defined FDR.
func calculateBenjamaniHochbergCutoff(contactScoreCache [][]int, s Settings, comparisonCountCache []int) []int {
	var cutoffCache = make([]int, len(contactScoreCache))
	var err error
	var currScore int
	var currQValue, currR, currP float64
	var out *fileio.EasyWriter
	var currRank int
	for i := range cutoffCache {
		cutoffCache[i] = s.MinCutoff
	}

	if s.FitStatsFile != "" {
		out = fileio.EasyCreate(s.FitStatsFile)
		_, err = fmt.Fprintf(out, "BinDistance\tR\tP\tCutoff\n")
		exception.PanicOnErr(err)
	}

	for currBinDistance := s.MinBinDistance; currBinDistance < len(contactScoreCache); currBinDistance++ {
		currRank = 0
		currR, currP = fit.ZeroTruncatedNegativeBinomial(contactScoreCache[currBinDistance], s.RStart, s.PStart, s.RStep, s.PStep)
		//fmt.Printf("Distance: %v. CurrR: %v. CurrP: %v.\n", currBinDistance, currR, currP)
		for currScore = len(contactScoreCache[currBinDistance]) - 1; currScore > s.MinCutoff; currScore-- {
			//fmt.Printf("CurrRank: %v. comparisonCountCache[currBinDistance]: %v.\n", currRank, comparisonCountCache[currBinDistance])
			currRank += contactScoreCache[currBinDistance][currScore]
			currQValue = (1 - numbers.NegativeBinomialCdf(float64(currScore), currR, currP)) * float64(comparisonCountCache[currBinDistance]) / float64(currRank)
			if !math.IsNaN(currQValue) && !math.IsInf(currQValue, 1) && !math.IsInf(currQValue, -1) {
				//fmt.Printf("CurrBinDistance: %v. CurrScore: %v. CurrQValue: %v.\n(1 - numbers.NegativeBinomialCdf(float64(currScore), currR, currP)): %v. float64(comparisonCountCache[currBinDistance]) / float64(currRank): %v. \n", currBinDistance, currScore, currQValue, 1-numbers.NegativeBinomialCdf(float64(currScore), currR, currP), float64(comparisonCountCache[currBinDistance])/float64(currRank))
			}
			if !math.IsNaN(currQValue) && !math.IsInf(currQValue, 1) && !math.IsInf(currQValue, -1) && currQValue > s.Fdr {
				cutoffCache[currBinDistance] = currScore
				break
			}
		}
		if s.FitStatsFile != "" {
			_, err = fmt.Fprintf(out, "%v\t%v\t%v\t%v\n", currBinDistance, currR, currP, cutoffCache[currBinDistance])
			exception.PanicOnErr(err)
		}
		if cutoffCache[currBinDistance] == s.MinCutoff {
			break
		}
	}

	if s.FitStatsFile != "" {
		err = out.Close()
		exception.PanicOnErr(err)
	}

	return cutoffCache
}

func usage() {
	fmt.Print(
		"strawToBedpe - Convert HiC contact maps in straw format to significant bedpe contact peak calls.\n" +
			"This program fits the distribution of contact scores for each distance to a zero-truncated negative\n" +
			"binomial distribution by regression evaluated by coordinate ascent. Then, significant peaks are reported\n" +
			"as significant outliers from these null contact distributions based on the total number of comparisons" +
			"and a user-defined false discovery rate.\n" +
			"The user must provide a txt file with the following tab separated fields:\n" +
			"file1.straw\tchrom\n" +
			"In the line above, the first field specifies a straw file, and the second field specifies\n" +
			"the chromosome name.\n" +
			"Inter-chromosomal straw files (where aChrom != bChrom) are not supported.\n" +
			"\n" +
			"Usage:\n" +
			"strawToBedpe [options] fileListAndChromNames.txt out.bedpe\n" +
			"options:\n")
	flag.PrintDefaults()
}

func main() {
	var expectedNumArgs int = 2
	var binSize *int = flag.Int("binSize", 5000, "Specify the binSize of input straw data.")
	var rStart *float64 = flag.Float64("rStart", 1.0, "Specify the starting parameter of R for coordinate ascent.")
	var pStart *float64 = flag.Float64("pStart", 0.5, "Specify the starting parameter of P for coordinate ascent.")
	var rStep *float64 = flag.Float64("rStep", 0.001, "Specify the step size for R for coordinate ascent.")
	var pStep *float64 = flag.Float64("pStep", 0.001, "Specify the step size for P for coordinate ascent.")
	var Fdr *float64 = flag.Float64("fdr", 0.05, "Set the false discovery rate for output peaks.")
	var minCutoff *int = flag.Int("minCutoff", 10, "Set the minimum number of reads required to call a significant peak.")
	var fitStatsFile *string = flag.String("fitStatsFile", "", "Write statistics about distribution fitting to an output text file.")
	var minBinDistance *int = flag.Int("minBinDistance", 0, "Specify the minimum bin distance required to call a significant peak.")

	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime | log.Lshortfile)
	flag.Parse()

	if len(flag.Args()) != expectedNumArgs {
		flag.Usage()
		log.Fatalf("Error: expecting %d arguments, but got %d\n",
			expectedNumArgs, len(flag.Args()))
	}

	fileList := flag.Arg(0)
	outFile := flag.Arg(1)

	s := Settings{
		FileList:       fileList,
		OutFile:        outFile,
		BinSize:        *binSize,
		RStart:         *rStart,
		PStart:         *pStart,
		RStep:          *rStep,
		PStep:          *pStep,
		MinCutoff:      *minCutoff,
		MinBinDistance: *minBinDistance,
		Fdr:            *Fdr,
		FitStatsFile:   *fitStatsFile,
	}

	strawToBedpe(s)
}
