// Command Group: "Data Conversion"

package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/convert"
	"github.com/vertgenlab/gonomics/fasta"
	"log"
)

func bedToFasta(fastaFile string, bedFile string, outfile string) {
	var records []bed.Bed = bed.Read(bedFile)
	var reference []fasta.Fasta = fasta.Read(fastaFile)
	var outlist []fasta.Fasta = convert.BedToFasta(records, reference)
	fasta.Write(outfile, outlist)
}

func usage() {
	fmt.Print(
		"bedToFasta - Extracts sequences from a fasta file from regions specified by an input bed.\n" +
			"Usage:\n" +
			"bedToFasta reference.fa intput.bed output.fa\n" +
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

	reference := flag.Arg(0)
	infile := flag.Arg(1)
	outfile := flag.Arg(2)

	bedToFasta(reference, infile, outfile)
}
