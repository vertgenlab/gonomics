package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/bed/bedpe"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/interval"
	"log"
)

func intervalContacts(bedpeFile string, inFile string, contactOutFile string) {
	var inIntervals = make([]interval.Interval, 0)
	var currOverlaps []interval.Interval
	var err error

	out := fileio.EasyCreate(contactOutFile)

	contactRecords := bedpe.Read(bedpeFile)
	inChan := interval.GoReadToChan(inFile)

	for i := range inChan {
		inIntervals = append(inIntervals, i)
	}
	inTree := interval.BuildTree(inIntervals)

	for _, i := range contactRecords {
		currOverlaps = interval.Query(inTree, i.A, "any")
		if len(currOverlaps) > 0 {
			bed.WriteBed(out, i.B) // if anything in the input set overlaps i.A, write i.B to the output.
		}
		currOverlaps = interval.Query(inTree, i.B, "any")
		if len(currOverlaps) > 0 {
			bed.WriteBed(out, i.A) // and vice versa
		}
	}

	err = out.Close()
	exception.PanicOnErr(err)
}

func usage() {
	fmt.Print("intervalContacts - Returns all regions that contact input genomic regions. " +
		"Contacts are specified by an input bedpe file.\n" +
		"Usage:\n" +
		"	intervalContacts [options] contacts.bedpe in.interval out.bed")
	flag.PrintDefaults()
}

func main() {
	var expectedNumArgs int = 3
	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime)

	if len(flag.Args()) != expectedNumArgs {
		flag.Usage()
		log.Fatalf("Error: expecting %d arguments, but got %d\n\n", expectedNumArgs, len(flag.Args()))
	}

	bedpeFile := flag.Arg(0)
	inFile := flag.Arg(1)
	contactOutFile := flag.Arg(2)

	intervalContacts(bedpeFile, inFile, contactOutFile)
}