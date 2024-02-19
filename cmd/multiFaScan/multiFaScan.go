// Command Group: "FASTA and Multi-FASTA Tools"

// Scan multiple Fasta alignment for a user-specified pattern (N for now) and report bed regions in reference sequence coordinates
package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/fasta"
	"log"
)

func multiFaScan(s Settings) {
	aln := fasta.Read(s.InFile)
	answerBedPos := fasta.ScanN(aln, s.QueryName)
	var currentBed bed.Bed
	var ans []bed.Bed
	for i := 0; i < len(answerBedPos); i++ {
		currentBed = bed.Bed{Chrom: s.Chrom, ChromStart: answerBedPos[i][0], ChromEnd: answerBedPos[i][1], Name: aln[0].Name, FieldsInitialized: 4} // Name field is reference species name
		ans = append(ans, currentBed)
	}
	bed.Write(s.OutFile, ans)
}

func usage() {
	fmt.Print(
		"multiFaScan - Scan multiple Fasta alignment for a user-specified pattern (N for now) and report bed regions in reference sequence coordinates.\n" +
			"Note that the reference sequence is assumed to be the first sequence in the multiple Fasta alignment.\n" +
			"Usage:\n" +
			"multiFaScan multi.fa out.fa\n" +
			"options:\n")
	flag.PrintDefaults()
}

type Settings struct {
	InFile    string
	OutFile   string
	QueryName string
	Chrom     string
}

func main() {
	var expectedNumArgs int = 2
	var queryName *string = flag.String("queryName", "", "Specify the name of the sequence to scan for the user-specified pattern (N for now). Set to the reference sequence by default")
	var chrom *string = flag.String("chrom", "chrom", "Specify the chromosome name of the reference species in the output bed file. Set to 'chrom' by default")
	var s Settings

	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime | log.Lshortfile)
	flag.Parse()

	if len(flag.Args()) != expectedNumArgs {
		flag.Usage()
		log.Fatalf("Error: expecting %d arguments, but got %d\n",
			expectedNumArgs, len(flag.Args()))
	}

	s = Settings{
		InFile:    flag.Arg(0),
		OutFile:   flag.Arg(1),
		QueryName: *queryName,
		Chrom:     *chrom,
	}

	multiFaScan(s)
}
