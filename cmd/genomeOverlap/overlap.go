package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/axt"
	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/interval"
	"log"
	"path"
	"strconv"
	"strings"
)

func usage() {
	fmt.Print(
		"genomeOverlap - Selects records in 'target' file that overlap with records in 'query' file.\n" +
			"Valid file types for target and query: Axt, Bed, Vcf, Sam \n" +
			"Usage:\n" +
			" genomeOverlap [options] \n" +
			"\t-query input.file\n" +
			"\t-target input.file \n" +
			"\t-out output.file \n\n" +
			"options:\n")
	flag.PrintDefaults()
}

func parseSelectRange(input string) *bed.Bed {
	words := strings.Split(input, ":")
	if len(words) != 2 {
		log.Fatalln("ERROR: Incorrect format in selectRange. Must be chr:start-end")
	}
	val := strings.Split(words[1], "-")
	if len(val) != 2 {
		log.Fatalln("ERROR: Incorrect format in selectRange. Must be chr:start-end")
	}
	start, _ := strconv.ParseInt(val[0], 10, 64)
	end, _ := strconv.ParseInt(val[1], 10, 64)
	answer := &bed.Bed{
		Chrom: words[0],
		ChromStart: start,
		ChromEnd: end}
	return answer
}

func goReadToIntervalChan(inputFile string) chan interval.Interval {
	answer := make(chan interval.Interval)
	go readToIntervalChan(inputFile, answer)
	return answer
}

func readToIntervalChan(inputFile string, send chan interval.Interval) {
	// How the file is read is dependent on the file extension
	filetype := path.Ext(inputFile)

	if filetype == ".gz" {
		// If terminal extension is ".gz" then trim off the gz and get the next extension
		filetype = path.Ext(inputFile[0 : len(inputFile)-len(filetype)])
	}

	switch filetype {
	case ".bed":
		receive := bed.GoReadToChan(inputFile)
		for val := range receive {
			send <- val
		}

	case ".axt":
		receive := axt.ReadToChan
	}
}

func genomeOverlap(query string, selectRange string, target string, out string, relationship string) {

	var inputQuery interval.Interval

	if selectRange != "" {
		inputQuery = parseSelectRange(selectRange)
	}
}

func main() {
	var helpRelationships *bool = flag.Bool("printRelationships", false, "Show a diagram of the valid interval relationships that can be tested for")
	var relationship *string = flag.String("relationship", "any", "Choose a specific relationships that target " +
		"and query records must fulfill to be reported. Use --printRelationships for more information.")
	var query *string = flag.String("query", "", "File containing query records. File will be read based on its file extension. See valid file types above.")
	var target *string = flag.String("target", "", "File containing target records that will be reported in output if they overlap a query record" +
		"File will be read based on its file extension. See valid file types above.")
	var out *string = flag.String("out", "", "Filename containing records from target that overlap query. File type will be the same as target file type")
	var selectRange *string = flag.String("selectRange", "", "Input a range and select all regions in target file that fall in this range. " +
		"Cannot be used at the same time as query. \n Format:\tchr:start-end\n NOTE: start and end must be 0 based (same as Bed format)")

	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime | log.Lshortfile)
	flag.Parse()

	if *helpRelationships {
		printRelationships()
		return
	}

	if !interval.TestValidRelationship(*relationship) {
		printRelationships()
		log.Fatalln("ERROR: Invalid relationship", *relationship)
	}

	if *query != "" && *selectRange != "" {
		log.Fatalln("ERROR: Cannot use both -query and -selectRange")
	}

	genomeOverlap(*query, *selectRange, *target, *out, *relationship)
}
