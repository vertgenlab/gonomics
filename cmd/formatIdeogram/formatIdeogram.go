// Command Group: "BED Tools"

// Generates ideogram txt file from a bed
package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/numbers/parse"
	"log"
	"strings"

	"github.com/vertgenlab/gonomics/fileio"
)

type IdeogramPoint struct {
	Chrom    string
	Position int
	Score    int
}

func formatIdeogram(inBed string, outTxt string, noScore bool) {
	var line string
	var doneReading bool = false
	var startNum, endNum, midpoint int
	var chrom string
	var outIdeogram []IdeogramPoint
	var err error

	file := fileio.EasyOpen(inBed)
	for line, doneReading = fileio.EasyNextRealLine(file); !doneReading; line, doneReading = fileio.EasyNextRealLine(file) {
		words := strings.Split(line, "\t")
		chrom = words[0]
		startNum = parse.StringToInt(words[1])
		endNum = parse.StringToInt(words[2])
		midpoint = (startNum + endNum) / 2

		outIdeogram = append(outIdeogram, IdeogramPoint{Chrom: chrom, Position: midpoint - 1, Score: 1})
		if noScore {
			outIdeogram = append(outIdeogram, IdeogramPoint{Chrom: chrom, Position: midpoint, Score: 10})
		} else {
			outIdeogram = append(outIdeogram, IdeogramPoint{Chrom: chrom, Position: midpoint, Score: parse.StringToInt(words[4])})
		}
		outIdeogram = append(outIdeogram, IdeogramPoint{Chrom: chrom, Position: midpoint + 1, Score: 1})
	}
	err = file.Close()
	exception.PanicOnErr(err)

	outfile := fileio.EasyCreate(outTxt)
	for i := 0; i < len(outIdeogram); i++ {
		_, err = fmt.Fprintf(outfile, "%s\t%v\t%v\n", outIdeogram[i].Chrom, outIdeogram[i].Position, outIdeogram[i].Score)
		exception.PanicOnErr(err)
	}
	err = outfile.Close()
	exception.PanicOnErr(err)
}

func usage() {
	fmt.Print(
		"formatIdeogram - Generates ideogram txt file from a bed.\n" +
			"Usage:\n" +
			"formatIdeogram input.bed output.txt\n" +
			"\n" +
			"How to use the output file of this program:\n" +
			"Run the program for a given bed file following the usage in the program itself.\n" +
			"To visualize the output file, go to the UCSC Genome Browser and navigate to the \n" +
			"\"Genome Graphs\" utility in the \"Tools\" tab.\n" +
			"You will then want to click on \"Upload\" and then add the output text file\n" +
			"generated by this program.\n" +
			"From there, you can use the GUI options on the Genome Graphs tool to produce\n" +
			"the desired visualization.\n" +
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
