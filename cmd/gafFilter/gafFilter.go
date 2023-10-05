// Command Group: "Ontology Tools"

// Filter gaf file.
package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/ontology/gaf"
	"log"
	"strings"
)

type Settings struct {
	InFile    string
	OutFile   string
	RemoveNot bool
}

func gafFilter(s Settings) {
	var pass bool

	gafChan, header := gaf.GoReadToChan(s.InFile)
	out := fileio.EasyCreate(s.OutFile)

	gaf.WriteHeaderToFileHandle(out, header)

	for curr := range gafChan {
		pass = true
		if s.RemoveNot && strings.Contains(curr.Qualifier, "NOT") {
			pass = false
		}

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
		"gafFilter\n" +
			"Usage:\n" +
			"gafFilter input.gaf output.gaf\n" +
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
