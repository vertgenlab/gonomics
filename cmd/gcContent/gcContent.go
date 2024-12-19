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

func gcContent(bedFile string, faFile string, outFile string, multiFaMode bool, species string) {

	var currRecordSeq []dna.Base
	var found bool
	var gc float64
	var annotated bed.Bed
	var speciesStart, speciesEnd int

	regionsChan := bed.GoReadToChan(bedFile)
	records := fasta.Read(faFile)
	recordsMap := fasta.ToMap(records)
	out := fileio.EasyCreate(outFile)

	for curr := range regionsChan {

		if multiFaMode {
			// when multiFa, each fasta sequence is a species, e.g. >human, >hca, assuming that the user specified the correct chromosome, e.g. this multi-fasta is chr1.fa
			currRecordSeq, found = recordsMap[species]
			if found {
				annotated = bed.Bed{Chrom: curr.Chrom, ChromStart: curr.ChromStart, ChromEnd: curr.ChromEnd, FieldsInitialized: 4}
				// need to convert bed region (reference sequence e.g. human RefPos), to AlnPos, then use AlnPos to index requested sequence (e.g. hca)
				speciesStart = fasta.RefPosToAlnPos(records[0], curr.ChromStart)
				speciesEnd = fasta.RefPosToAlnPos(records[0], curr.ChromEnd)
				gc = dna.GCContent(currRecordSeq[speciesStart:speciesEnd])
				annotated.Name = fmt.Sprintf("%e", gc)
				bed.WriteBed(out, annotated)
			} else {
				log.Fatalf("Error: multiFaMode. Requested species (%s) was not found as a fasta record name in the input multi-fasta file\n", species)
			}

		} else {
			// when not a multiFa, assuming genome.fa, where each fasta sequence is a chromosome, e.g. >chr1
			currRecordSeq, found = recordsMap[curr.Chrom]
			if found {
				annotated = bed.Bed{Chrom: curr.Chrom, ChromStart: curr.ChromStart, ChromEnd: curr.ChromEnd, FieldsInitialized: 4}
				gc = dna.GCContent(currRecordSeq[curr.ChromStart:curr.ChromEnd])
				annotated.Name = fmt.Sprintf("%e", gc)
				bed.WriteBed(out, annotated)
			} else {
				log.Fatalf("Error: bed region chrom (%s) was not found as a fasta record name in the input fasta file\n", curr.Chrom)
			}
		}

	}

	err := out.Close()
	exception.PanicOnErr(err)
}

func usage() {
	fmt.Print(
		"gcContent - Calculates the gc content of fasta sequences for all regions specified in a bed file.\n" +
			"Outputs a new bed file with 4 columns: Chrom, ChromStart, ChromEnd, gcContent.\n" +
			"Usage:\n" +
			"gcContent in.bed in.fa mult.fa out.bed\n" +
			"options:\n")
	flag.PrintDefaults()
}

func main() {
	var expectedNumArgs int = 3
	var multiFaMode *bool = flag.Bool("multiFaMode", false, "When in multiFa mode, each fasta sequence is a species, e.g. >human, >hca, assuming that the user specified the correct chromosome, e.g. this multi-fasta is chr1.fa. When not in multiFa mode, assuming genome.fa, where each fasta sequence is a chromosome, e.g. >chr1.")
	var species *string = flag.String("multiFaSpecies", "", "When in multiFa mode, must specify the species to calculate gc content for.")

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

	gcContent(bedFile, faFile, outFile, *multiFaMode, *species)
}
