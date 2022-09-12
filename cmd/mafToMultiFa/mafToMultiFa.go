// Command Group: "Data Conversion"

package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/maf"
	"log"
)

func mafToMultiFa(inMaf string, inFa string, inSpeciesList string, outMultiFa string, noMask bool) {
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

	fasta.Write(outMultiFa, aln)
}

func usage() {
	fmt.Print(
		"mafToMultiFa - convert a maf alignment into a multiple fasta alignment\n" +
			"Usage:\n" +
			" mafToMultiFa input.maf reference.fa species.list output.fa\n" +
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
	outMultiFaFile := flag.Arg(3)

	mafToMultiFa(mafFile, refFile, speciesFile, outMultiFaFile, *noMask)
}
