package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/bed"
	"log"
	"strings"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/common"
)

func bedFilter(infile string, outfile string, minScore int64, maxScore int64, minLength int64, maxLength int64, minStart int64, maxStart int64, minEnd int64, maxEnd int64, chrom string) {
	var outlist []*bed.Bed
	var line string
	var startNum, endNum, length int64
	var doneReading bool = false
	var pass bool = false
	var current *bed.Bed
	file := fileio.EasyOpen(infile)
	defer file.Close()

	for line, doneReading = fileio.EasyNextRealLine(file); !doneReading; line, doneReading = fileio.EasyNextRealLine(file) {
		pass = true
		words := strings.Split(line, "\t")
		startNum = common.StringToInt64(words[1])
		endNum = common.StringToInt64(words[2])
		length = endNum - startNum
		if len(words) > 4{
			if common.StringToInt64(words[4]) < minScore {
				pass = false
		} 
			if maxScore != 0 {
				if common.StringToInt64(words[4]) > maxScore {
					pass = false
				}
			} 

		}
		
		if length < minLength {
			pass = false
		} 
		if maxLength != 0 {
			if length > maxLength {
				pass = false
			}
		} 
		if startNum < minStart {
			pass = false
		}
		if maxStart != 0 {
			if startNum > maxStart {
				pass = false
			}
		} 
		if endNum < minEnd {
			pass = false
		}
		if maxStart != 0 {
			if endNum > maxEnd {
				pass = false
			}
		} 
		if chrom != "" {
			if words[0] != chrom {
				pass = false
			}
		}

		if pass {
			current = &bed.Bed{Chrom: words[0], ChromStart: startNum, ChromEnd: endNum, Name: words[1], Score: common.StringToInt64(words[4])}
			outlist = append(outlist, current)
		}
	}
	bed.Write(outfile, outlist, 5)
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
	var minScore *int64 = flag.Int64("minScore", 0, "Specifies the minimum score in the fourth field.")
	var maxScore *int64 = flag.Int64("maxScore", 0, "Specifies the maximum score in the fourth field.")
	var minLength *int64 = flag.Int64("minLength", 0, "Specifies the minimum length of the region.")
	var maxLength *int64 = flag.Int64("maxLength", 0, "Specifies the maximum length of the region.")
	var minStart *int64 = flag.Int64("minStart", 0, "Specifies the minimum starting position of the region.")
	var maxStart *int64 = flag.Int64("maxStart", 0, "Specifies the maximum starting position of the region.")
	var minEnd *int64 = flag.Int64("minEnd", 0, "Specifies the minimum ending position of the region.")
	var maxEnd *int64 = flag.Int64("maxEnd", 0, "Specifies the maximum ending position of the region.")
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
