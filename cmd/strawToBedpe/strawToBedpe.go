// Command Group: "Data Conversion"

// Convert HiC contact maps in straw format to bedpe contact peak calls
package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/hic"
	"github.com/vertgenlab/gonomics/numbers"
	"log"
	"math"
	"strings"
)

func strawToBedpe(fileList string, outFile string, binSize int) {
	//var thisRec bedpe.BedPe
	//var err error
	var currStrawChan <-chan hic.Straw
	var words []string
	//var currAChrom, currBChrom string
	var binDistance int
	var currDistance float64
	var contactScoreCache [][]int = make([][]int, 0)
	var searchSpaceMins = make(map[string]int)
	var searchSpaceMaxes = make(map[string]int)
	var tmpContactScoreCache [][]int
	var tmpColumn []int
	var firstTime, foundInMap bool

	lines := fileio.Read(fileList)

	for i := range lines {
		words = strings.Split(lines[i], "\t")
		currStrawChan = hic.GoReadToChan(words[0])
		firstTime = true
		for s := range currStrawChan {
			if firstTime {
				// if we've already had a bedpe file with information from this chromosome.
				if _, foundInMap = searchSpaceMins[words[1]]; foundInMap {
					searchSpaceMins[words[1]] = numbers.Min(searchSpaceMins[words[1]], numbers.Min(s.Bin1Start, s.Bin2Start))
					searchSpaceMaxes[words[1]] = numbers.Max(searchSpaceMaxes[words[1]], numbers.Max(s.Bin1Start, s.Bin2Start))
				} else {
					searchSpaceMins[words[1]] = numbers.Min(s.Bin1Start, s.Bin2Start)
					searchSpaceMaxes[words[1]] = numbers.Max(s.Bin1Start, s.Bin2Start)
				}
			}

			if s.Bin1Start < searchSpaceMins[words[1]] {
				searchSpaceMins[words[1]] = s.Bin1Start
			}
			if s.Bin2Start < searchSpaceMins[words[1]] {
				searchSpaceMins[words[1]] = s.Bin2Start
			}
			if s.Bin1Start > searchSpaceMaxes[words[1]] {
				searchSpaceMaxes[words[1]] = s.Bin1Start
			}
			if s.Bin2Start > searchSpaceMaxes[words[1]] {
				searchSpaceMaxes[words[1]] = s.Bin2Start
			}

			currDistance = math.Abs(float64(s.Bin1Start) - float64(s.Bin2Start))
			if int(currDistance)%binSize != 0 {
				log.Fatalf("Error: Distance between two straw ends: %v is not a multiple of the bin size: %v.\n", currDistance, binSize)
			}
			binDistance = int(currDistance) / binSize
			// if we need more room in the cache rows, extend the cache
			if binDistance+1 > len(contactScoreCache) {
				tmpContactScoreCache = make([][]int, binDistance+1)
				copy(tmpContactScoreCache, contactScoreCache)
				contactScoreCache = tmpContactScoreCache
			}

			if contactScoreCache[binDistance] == nil {
				contactScoreCache[binDistance] = make([]int, 1)
			}
			// if we need more room in the cache columns, extend the cache
			if s.ContactScore+1 > len(contactScoreCache[binDistance]) {
				tmpColumn = make([]int, s.ContactScore+1)
				copy(tmpColumn, contactScoreCache[binDistance])
				contactScoreCache[binDistance] = tmpColumn
			}
			contactScoreCache[binDistance][s.ContactScore]++
		}
	}

	var currChrom string
	var currNumWindows int
	var currSumNonZeroWindows int
	for currBinDistance := range contactScoreCache {
		// get the total numbers of windows at this binDistance by ranging across all chromosomes
		for currChrom = range searchSpaceMins {
			// the current number of windows on a chromosome is the size of the search space (the highest window minus the lowest window)
			// divided by the binSize (this is the number of bins on the chromosome), minus the binDistance ( as the right most
			// possible start positions will not form a valid window).
			currNumWindows = (searchSpaceMaxes[currChrom]-searchSpaceMins[currChrom])/binSize - currBinDistance
		}
		currSumNonZeroWindows = 0
		for j := range contactScoreCache[currBinDistance] {
			currSumNonZeroWindows += contactScoreCache[currBinDistance][j]
		}
		if contactScoreCache[currBinDistance] == nil {
			contactScoreCache[currBinDistance] = make([]int, 1)
		}
		contactScoreCache[currBinDistance][0] = currNumWindows - currSumNonZeroWindows
	}

	out := fileio.EasyCreate("5kb.txt")
	_, err := fmt.Fprintf(out, "Score\tFrequency\n")
	exception.PanicOnErr(err)

	for i := range contactScoreCache[1] {
		_, err = fmt.Fprintf(out, "%v\t%v\n", i, contactScoreCache[1][i])
		exception.PanicOnErr(err)
	}

	err = out.Close()
	exception.PanicOnErr(err)

	out = fileio.EasyCreate("0kb.txt")
	_, err = fmt.Fprintf(out, "Score\tFrequency\n")
	exception.PanicOnErr(err)

	for i := range contactScoreCache[0] {
		_, err = fmt.Fprintf(out, "%v\t%v\n", i, contactScoreCache[0][i])
		exception.PanicOnErr(err)
	}

	err = out.Close()
	exception.PanicOnErr(err)

	out = fileio.EasyCreate("10kb.txt")
	_, err = fmt.Fprintf(out, "Score\tFrequency\n")
	exception.PanicOnErr(err)

	for i := range contactScoreCache[2] {
		_, err = fmt.Fprintf(out, "%v\t%v\n", i, contactScoreCache[2][i])
		exception.PanicOnErr(err)
	}

	err = out.Close()
	exception.PanicOnErr(err)

	//fmt.Printf("ChromStarts: %v.\n", searchSpaceMins)
	//fmt.Printf("ChromEnds: %v.\n", searchSpaceMaxes)
	//fmt.Printf("ContactScoreCache:\n%v", contactScoreCache)

	/*

		straw := hic.GoReadToChan(strawFile)
		out := fileio.EasyCreate(outFile)

		if interChrom == "" {
			for s := range straw {
				thisRec = bedpe.BedPe{A: bed.Bed{Chrom: chrom, ChromStart: s.Bin1Start, ChromEnd: s.Bin1Start + binSize, Score: s.ContactScore, FieldsInitialized: 8}, B: bed.Bed{Chrom: chrom, ChromStart: s.Bin2Start, ChromEnd: s.Bin2Start + binSize, Score: s.ContactScore, FieldsInitialized: 8}}
				bedpe.WriteToFileHandle(out, thisRec)
			}
		} else {
			for s := range straw {
				thisRec = bedpe.BedPe{A: bed.Bed{Chrom: chrom, ChromStart: s.Bin1Start, ChromEnd: s.Bin1Start + binSize, Score: s.ContactScore, FieldsInitialized: 8}, B: bed.Bed{Chrom: interChrom, ChromStart: s.Bin2Start, ChromEnd: s.Bin2Start + binSize, Score: s.ContactScore, FieldsInitialized: 8}}
				bedpe.WriteToFileHandle(out, thisRec)
			}
		}

		err = out.Close()
		exception.PanicOnErr(err)
	*/
}

func usage() {
	fmt.Print(
		"strawToBedpe - Convert HiC contact maps in straw format to bedpe contact peak calls.\n" +
			"The user must provide a txt file with the following tab separated fields:\n" +
			"file1.straw\tchrom\n" +
			"In the line above, the first field specifies a straw file, and the second field specifies\n" +
			"the chromosome name.\n" +
			"Inter-chromosomal straw files (where aChrom != bChrom), are not currently supported.\n" +
			"\n" +
			"Usage:\n" +
			"strawToBedpe [options] fileListAndChromNames.txt out.bedpe\n" +
			"options:\n")
	flag.PrintDefaults()
}

func main() {
	var expectedNumArgs int = 2
	var binSize *int = flag.Int("binSize", 5000, "Specify the binSize of input straw data.")

	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime | log.Lshortfile)
	flag.Parse()

	if len(flag.Args()) != expectedNumArgs {
		flag.Usage()
		log.Fatalf("Error: expecting %d arguments, but got %d\n",
			expectedNumArgs, len(flag.Args()))
	}

	fileListAndChromFile := flag.Arg(0)
	outFile := flag.Arg(1)

	strawToBedpe(fileListAndChromFile, outFile, *binSize)
}
