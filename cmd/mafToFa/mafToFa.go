package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/maf"
	"log"
)

func mafToFa(inMaf string, inFa string, inSpeciesList string, outFa string, noMask bool) {
	mafRecords := maf.Read(inMaf)

	refFasta := fasta.Read(inFa)
	if len(refFasta) != 1 {
		log.Fatalf("Error: expecting input fasta to be a single record, but file has %d records\n", len(refFasta))
	}

	speciesList := fileio.Read(inSpeciesList)

	aln := maf.ToFasta(mafRecords, refFasta[0], speciesList)

	if noMask {
		fasta.AllToUpper(aln)
	}

	fasta.Write(outFa, aln)
}

func usage() {
	fmt.Print(
		"mafToFa - convert a maf alignment into a fasta alignment\n" +
			"Usage:\n" +
			" mafToFa input.maf reference.fa species.list output.fa\n" +
			"options:\n")
	flag.PrintDefaults()
}

func main() {
	var expectedNumArgs int = 4
	var noMask *bool = flag.Bool("noMask", false, "convert all bases to uppercase in the output")
	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime | log.Lshortfile)
	flag.Parse()

	if len(flag.Args()) != expectedNumArgs {
		flag.Usage()
		log.Fatalf("Error: expecting %d arguments, but got %d\n", expectedNumArgs, len(flag.Args()))
	}

	mafFile := flag.Arg(0)
	refFile := flag.Arg(1)
	speciesFile := flag.Arg(2)
	outFile := flag.Arg(3)

	mafToFa(mafFile, refFile, speciesFile, outFile, *noMask)
}
