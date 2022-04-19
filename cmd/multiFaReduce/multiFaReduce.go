// Command Group: "FASTA and Multi-FASTA Tools"

package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/fasta"
	"log"
)

func makeEmptyFastas(records []fasta.Fasta) []fasta.Fasta {
	var answer []fasta.Fasta = make([]fasta.Fasta, len(records))
	for i := range records {
		answer[i].Name = records[i].Name
	}
	return answer
}

func colVariable(records []fasta.Fasta, colIdx int) bool {
	var i int
	var firstBase dna.Base

	firstBase = records[0].Seq[colIdx]
	for i = 1; i < len(records); i++ {
		if records[i].Seq[colIdx] != firstBase {
			return true
		}
	}
	return false
}

func addCopyOfCol(from []fasta.Fasta, alnIdx int, to []fasta.Fasta) {
	for i := range from {
		to[i].Seq = append(to[i].Seq, from[i].Seq[alnIdx])
	}
}

func mfaReduce(inFilename, outFilename string) {
	records := fasta.Read(inFilename)

	answer := makeEmptyFastas(records)

	for i := range records[0].Seq {
		if colVariable(records, i) {
			addCopyOfCol(records, i, answer)
		}
	}

	fasta.Write(outFilename, answer)
}

func usage() {
	fmt.Print(
		"mfaReduce - mfaReduce removes all columns in an multi fasta alignment that are not variable\n" +
			"Usage:\n" +
			" mfaReduce input.mfa output.mfa\n" +
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
		log.Fatalf("Error: expecting %d arguments, but got %d\n", expectedNumArgs, len(flag.Args()))
	}

	inFile := flag.Arg(0)
	outFile := flag.Arg(1)

	mfaReduce(inFile, outFile)
}
