// Command Group:

// Output a subset of a gaf file without "NOT" entries in qualifier field
package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/ontology/gaf"
	"log"
)

type Settings struct {
	InFile    string
	OutFile   string
	RemoveNot bool
}

func gafFilter(s Settings) {
	var pass bool = false

	gafChan, _ := gaf.GoReadToChan(s.InFile)
	out := fileio.EasyCreate(s.OutFile)

	for curr := range gafChan {
		pass = true
		/*if s.Chrom != "" {
			if curr.Chrom != s.Chrom {
				pass = false
			}
		}*/

		fmt.Printf("curr.Qualifier: %s", curr.Qualifier)

		if pass {
			gaf.WriteGaf(out, curr)
		}
	}
	var err error
	err = out.Close()
	exception.PanicOnErr(err)
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
	var removeNot *bool = flag.Bool("removeNot", false, "Specifies remove entries where qualifier contains NOT.")

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
	s := Settings{
		InFile:    infile,
		OutFile:   outfile,
		RemoveNot: *removeNot,
	}
	gafFilter(s)
}
