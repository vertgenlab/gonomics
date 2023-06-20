package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/fileio"
	"log"
)

func functionForTheTool(b string, f string, output string) {
	bedSeq := make([][]dna.Base, len(b))
	aaSeq := make([][]dna.AminoAcid, len(b))
	var outRecord []string
	bedRecords := bed.Read(b)
	fastaRecord := fasta.Read(f)

	for i := range bedRecords {
		for j := range fastaRecord[0].Seq { //this assumes one fasta record in the file
			if bedRecords[i].ChromStart <= j && j < bedRecords[i].ChromEnd {
				bedSeq[i] = append(bedSeq[i], fastaRecord[0].Seq[j])
			}
		}
		aaSeq[i] = append(aaSeq[i], dna.TranslateSeq(bedSeq[i])...)
	}

	for k := range aaSeq {
		for l := range aaSeq[k] {
			outRecord = append(outRecord, dna.AminoAcidToString(aaSeq[k][l]))
		}
	}

	fileio.Write(output, outRecord)

}

func usage() {
	fmt.Print(
		"newTool - takes a bed and fasta and converts the bed sequences into amino acide sequences\n" +
			"Usage:\n" +
			"newTool bedFile fastaFile\n" +
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
	fastaFile := flag.Arg(1)
	outFile := flag.Arg(2)

	functionForTheTool(bedFile, fastaFile, outFile)
}
