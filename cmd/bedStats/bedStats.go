//Command Group: BED Tools

package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/bed/bedpe"
	"log"
)

type Settings struct {
	InContacts string
	InTestBed  string
	OutMatched string
}

func bedStats(s Settings) {
	var freq float64
	var matchedBeds []bed.Bed
	trueContacts := bedpe.Read(s.InContacts)
	testBeds := bed.Read(s.InTestBed)
	freq, matchedBeds = bedpe.GeneAssignmentCheck(trueContacts, testBeds)

	bed.Write(s.OutMatched, matchedBeds)
	fmt.Println(freq)
}

func usage() {
	fmt.Print(
		"bedStats - takes a bedpe containing tru contacts from empirical data which assigns regions of the genome \n" +
			"to putative target genes and compares an output from assignGenomeSpace command to determine how \n" +
			"accurate the nearest gene in 3d space calculation was. Writes out a frequency of correct assignments \n" +
			"and outputs a bed containing the matching regions.\n" +
			"Usage:\n" +
			"bedStats true.bedpe test.bed matched.bed\n" +
			"options:\n")
	flag.PrintDefaults()
}
func main() {
	var expectedNumArgs int = 3

	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime | log.Lshortfile)
	flag.Parse()
	if len(flag.Args()) != expectedNumArgs {
		flag.Usage()
		log.Fatalf("Error: expecting %d arguments, but got %d\n",
			expectedNumArgs, len(flag.Args()))
	}
	inContacts := flag.Arg(0)
	inTest := flag.Arg(1)
	matchedBeds := flag.Arg(2)

	s := Settings{
		InContacts: inContacts,
		InTestBed:  inTest,
		OutMatched: matchedBeds,
	}
	bedStats(s)
}
