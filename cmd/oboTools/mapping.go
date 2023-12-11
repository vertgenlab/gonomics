package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/ontology/obo"
	"log"
	"os"
)

type MappingSettings struct {
	InFile  string
	OutFile string
	Force   bool
}

func mappingUsage(mappingFlags *flag.FlagSet) {
	fmt.Print("oboTools mapping - a subcommand to create tabular summaries of Obo ontologies.\n" +
		"Usage:\n" +
		"oboTools mapping in.obo out.txt\n" +
		"Options:\n")
	mappingFlags.PrintDefaults()
}

func parseMappingArgs() {
	var err error
	var expectedNumArgs int = 2 //excludes the subcommand

	mappingFlags := flag.NewFlagSet("mapping", flag.ExitOnError)
	force := mappingFlags.Bool("force", false, "(Advanced) Force the program to ignore parsing warnings, such as duplicate fields or missing parent nodes.")
	err = mappingFlags.Parse(os.Args[2:])
	mappingFlags.Usage = func() { mappingUsage(mappingFlags) }
	exception.PanicOnErr(err)

	if len(mappingFlags.Args()) != expectedNumArgs {
		mappingFlags.Usage()
		log.Fatalf("Error: Expected %v arguments. Found: %v.\n", expectedNumArgs, len(mappingFlags.Args()))
	}

	inFile := mappingFlags.Arg(0)
	outFile := mappingFlags.Arg(1)
	s := MappingSettings{
		InFile:  inFile,
		OutFile: outFile,
		Force:   *force,
	}

	OboToolsMapping(s)
}

func OboToolsMapping(s MappingSettings) {
	var err error
	records, _ := obo.Read(s.InFile, s.Force)
	out := fileio.EasyCreate(s.OutFile)
	for i := range records {
		_, err = fmt.Fprintf(out, "%s\t%s\n", records[i].Id, records[i].Name)
		exception.PanicOnErr(err)
	}

	err = out.Close()
	exception.PanicOnErr(err)
}
