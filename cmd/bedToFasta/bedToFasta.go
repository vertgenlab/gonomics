// Command Group: "Data Conversion"

package main

import (
	"flag"
	"fmt"
	"log"

	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/convert"
	"github.com/vertgenlab/gonomics/fasta"
)

func bedToFasta(fastaFile string, bedFile string, outfile string, revComp bool) {
	var records []bed.Bed = bed.Read(bedFile)
	var reference []fasta.Fasta = fasta.Read(fastaFile)
	var fastaEntry fasta.Fasta
	var outlist []fasta.Fasta

	for i := range records {
		fastaEntry = convert.SingleBedToFasta(records[i], reference)
		if (revComp == true) && (records[i].Strand == bed.Negative) {
			fasta.ReverseComplement(fastaEntry)
			fastaEntry.Name = fastaEntry.Name + "_RevComp"
		}
		outlist = append(outlist, fastaEntry)
	}
	fasta.Write(outfile, outlist)
}

func usage() {
	fmt.Print(
		"bedToFasta - Extracts sequences from a fasta file from regions specified by an input bed.\n" +
			"Usage:\n" +
			"bedToFasta reference.fa input.bed output.fa\n" +
			"options:\n")
	flag.PrintDefaults()
}

func main() {
	var expectedNumArgs int = 3
	var revComp *bool = flag.Bool("revComp", false, "Reverse complement fasta output if the bed strand is negative. Will append '_RevComp' to fasta name.")

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

	bedToFasta(reference, infile, outfile, *revComp)
}
