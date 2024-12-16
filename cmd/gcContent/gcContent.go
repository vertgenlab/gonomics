// Command Group: "FASTA and Multi-FASTA Tools"

// Calculates the gc content of fasta sequences for all regions specified in a bed file
package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/fileio"
	"log"
)

func gcContent(bedFile string, faFile string, outFile string) {

	var currRecordSeq []dna.Base
	var found bool
	var gc float64

	regionsChan := bed.GoReadToChan(bedFile)
	records := fasta.Read(faFile)
	recordsMap := fasta.ToMap(records)
	out := fileio.EasyCreate(outFile)

	for curr := range regionsChan {

		currRecordSeq, found = recordsMap[curr.Chrom]

		if found {
			gc = dna.GCContent(currRecordSeq[curr.ChromStart:curr.ChromEnd])
			curr.Annotation = append(curr.Annotation, fmt.Sprintf("%e", gc))
			bed.WriteBed(out, curr)
		} else {
			log.Fatalf("Error: bed region chrom (%s) was not found as a fasta record name in the input fasta file\n", curr.Chrom)
		}

	}

	err := out.Close()
	exception.PanicOnErr(err)
}

func usage() {
	fmt.Print(
		"gcContent - Calculates the gc content of fasta sequences for all regions specified in a bed file.\n" +
			"Usage:\n" +
			"gcContent in.bed in.fa mult.fa out.bed\n" +
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

	bedFile := flag.Arg(0)
	faFile := flag.Arg(1)
	outFile := flag.Arg(2)

	gcContent(bedFile, faFile, outFile)
}
