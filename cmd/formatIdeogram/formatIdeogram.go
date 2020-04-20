package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/common"
	"github.com/vertgenlab/gonomics/fileio"
	"log"
	"strings"
)

type IdeogramPoint struct {
	Chrom    string
	Position int64
	Score    int64
}

func formatIdeogram(inBed string, outTxt string, noScore bool) {
	var line string
	var doneReading bool = false
	var startNum, endNum, midpoint int64
	var chrom string
	var outIdeogram []*IdeogramPoint
	var err error

	file := fileio.EasyOpen(inBed)
	defer file.Close()

	for line, doneReading = fileio.EasyNextRealLine(file); !doneReading; line, doneReading = fileio.EasyNextRealLine(file) {
		words := strings.Split(line, "\t")
		chrom = words[0]
		startNum = common.StringToInt64(words[1])
		endNum = common.StringToInt64(words[2])
		midpoint = (startNum + endNum) / int64(2)

		outIdeogram = append(outIdeogram, &IdeogramPoint{Chrom: chrom, Position: midpoint - int64(1), Score: int64(1)})
		if noScore {
			outIdeogram = append(outIdeogram, &IdeogramPoint{Chrom: chrom, Position: midpoint, Score: 10})
		} else {
			outIdeogram = append(outIdeogram, &IdeogramPoint{Chrom: chrom, Position: midpoint, Score: common.StringToInt64(words[4])})
		}
		outIdeogram = append(outIdeogram, &IdeogramPoint{Chrom: chrom, Position: midpoint + int64(1), Score: int64(1)})
	}

	outfile := fileio.EasyCreate(outTxt)
	defer outfile.Close()

	for i := 0; i < len(outIdeogram); i++ {
		_, err = fmt.Fprintf(outfile, "%s\t%v\t%v\n", outIdeogram[i].Chrom, outIdeogram[i].Position, outIdeogram[i].Score)
		common.ExitIfError(err)
	}
}

func usage() {
	fmt.Print(
		"formatIdeogram - Generates ideogram txt file from a bed.\n" +
			"Usage:\n" +
			"formatIdeogram input.bed output.txt\n" +
			"options:\n")
	flag.PrintDefaults()
}

func main() {
	var expectedNumArgs int = 2
	var noScore *bool = flag.Bool("noScore", false, "Used for bed files without a score column.")
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

	formatIdeogram(infile, outfile, *noScore)
}
