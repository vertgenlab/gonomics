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
		"bedpeFilter\n" +
			"Usage:\n" +
			"bedpeFilter input.bedpe output.bedpe\n" +
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
