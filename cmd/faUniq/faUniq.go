package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/fasta"
	"log"
)

func faUniq(inFile string, outFile string) {
	allSeqs := fasta.Read(inFile)

	seqsSeen := make(map[string]int, 0) // This maps the DNA sequence converted to a string (key) to the index in uniqueSeqs (value)
	seqsSeen[dna.BasesToString(allSeqs[0].Seq)] = 0

	uniqueSeqs := make([]fasta.Fasta, 1)
	uniqueSeqs[0] = allSeqs[0]

	var currSeq string
	var seqExists bool
	var val int
	for i := 1; i < len(allSeqs); i++ {
		currSeq = dna.BasesToString(allSeqs[i].Seq)
		val, seqExists = seqsSeen[currSeq]
		if !seqExists {
			seqsSeen[currSeq] = len(uniqueSeqs)
			uniqueSeqs = append(uniqueSeqs, allSeqs[i])
		} else {
			uniqueSeqs[val].Name = uniqueSeqs[val].Name + "; " + allSeqs[i].Name
		}
	}

	fasta.Write(outFile, uniqueSeqs)
}

func usage() {
	fmt.Print(
		"faUniq - Pull unique sequences from a fasta file.\n" +
			"Usage:\n" +
			"faUniq inputFile.fa outputFile.fa\n" +
			"options:\n")
	flag.PrintDefaults()
}

func main() {
	var expectedNumArgs int = 2

	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime | log.Lshortfile)
	flag.Parse()

	if len(flag.Args()) != expectedNumArgs {
		flag.Usage()
		log.Fatalf("Error: expecting %d arguments, but got %d\n",
			expectedNumArgs, len(flag.Args()))
	}

	inputFile := flag.Arg(0)
	outputFile := flag.Arg(1)

	faUniq(inputFile, outputFile)
}
