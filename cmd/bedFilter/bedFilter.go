// Command Group: "BED Tools"

package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/common"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/numbers"
	"log"
	"strings"
)

func bedFilter(infile string, outfile string, minScore int, maxScore int, minLength int, maxLength int, minStart int, maxStart int, minEnd int, maxEnd int, chrom string) {
	var outlist []bed.Bed
	var line string
	var startNum, endNum, length int
	var doneReading bool = false
	var pass bool = false
	var numFields int
	var current bed.Bed
	file := fileio.EasyOpen(infile)
	defer file.Close()

	for line, doneReading = fileio.EasyNextRealLine(file); !doneReading; line, doneReading = fileio.EasyNextRealLine(file) {
		pass = true
		words := strings.Split(line, "\t")
		startNum = common.StringToInt(words[1])
		endNum = common.StringToInt(words[2])
		length = endNum - startNum
		if len(words) > 4 {
			if common.StringToInt(words[4]) < minScore {
				pass = false
			}
			if common.StringToInt(words[4]) > maxScore {
				pass = false
			}
		}
		if length < minLength {
			pass = false
		}
		if length > maxLength {
			pass = false
		}
		if startNum < minStart {
			pass = false
		}
		if startNum > maxStart {
			pass = false
		}
		if endNum < minEnd {
			pass = false
		}
		if endNum > maxEnd {
			pass = false
		}
		if chrom != "" {
			if words[0] != chrom {
				pass = false
			}
		}
		if pass {
			if len(words) == 3 {
				current = bed.Bed{Chrom: words[0], ChromStart: startNum, ChromEnd: endNum}
				numFields = numbers.Max(3, numFields)
			} else if len(words) == 4 {
				current = bed.Bed{Chrom: words[0], ChromStart: startNum, ChromEnd: endNum, Name: words[3]}
				numFields = numbers.Max(4, numFields)
			} else {
				current = bed.Bed{Chrom: words[0], ChromStart: startNum, ChromEnd: endNum, Name: words[3], Score: common.StringToInt(words[4])}
				numFields = numbers.Max(5, numFields)
			}
			outlist = append(outlist, current)
		}
	}
	bed.Write(outfile, outlist)
}

func usage() {
	fmt.Print(
		"bedFilter\n" +
			"Usage:\n" +
			"bedFilter input.bed output.bed\n" +
			"options:\n")
	flag.PrintDefaults()
}

func main() {
	var expectedNumArgs int = 2
	var minScore *int = flag.Int("minScore", 0, "Specifies the minimum score in the fourth field.")
	var maxScore *int = flag.Int("maxScore", numbers.MaxInt, "Specifies the maximum score in the fourth field.")
	var minLength *int = flag.Int("minLength", 0, "Specifies the minimum length of the region.")
	var maxLength *int = flag.Int("maxLength", numbers.MaxInt, "Specifies the maximum length of the region.")
	var minStart *int = flag.Int("minStart", 0, "Specifies the minimum starting position of the region.")
	var maxStart *int = flag.Int("maxStart", numbers.MaxInt, "Specifies the maximum starting position of the region.")
	var minEnd *int = flag.Int("minEnd", 0, "Specifies the minimum ending position of the region.")
	var maxEnd *int = flag.Int("maxEnd", numbers.MaxInt, "Specifies the maximum ending position of the region.")
	var chrom *string = flag.String("chrom", "", "Specifies the chromosome name.")

	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime | log.Lshortfile)
	flag.Parse()

	if len(flag.Args()) != expectedNumArgs {
		flag.Usage()
		log.Fatalf("Error: expecting %d arguments, but got %d\n",
			expectedNumArgs, len(flag.Args()))
	}

	infile := flag.Arg(0)
	outfile := flag.Arg(1)

	bedFilter(infile, outfile, *minScore, *maxScore, *minLength, *maxLength, *minStart, *maxStart, *minEnd, *maxEnd, *chrom)
}
