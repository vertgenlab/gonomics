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

// bedToAminoAcid is  complete function made from gonomics packages, it takes a bed and a fasta and returns the amino acid sequence of all the beds together.
func bedToAminoAcid(b string, f string, output string) {
	bedSeq := make([][]dna.Base, len(b))     //this creates a place to store all the bed regions' dna sequences
	aaSeq := make([][]dna.AminoAcid, len(b)) //this creates a place for all the amino acid sequences that correspond to the dna sequences
	var outRecord []string                   //this will be the variable that stores the output that will be written to the file
	bedRecords := bed.Read(b)                //reading the bed
	fastaRecord := fasta.Read(f)             //reading the fasta

	for i := range bedRecords { //for each bedRecord
		for j := range fastaRecord[0].Seq { //go through each base in the fasta record, assumes one fasta record in the file
			if bedRecords[i].ChromStart <= j && j < bedRecords[i].ChromEnd { //find locations where the coordinates in the fasta record match the bed coordinates
				bedSeq[i] = append(bedSeq[i], fastaRecord[0].Seq[j]) //add the sequence for the bed coordinates to the variable for the dna sequences
			}
		}
		aaSeq[i] = append(aaSeq[i], dna.TranslateSeq(bedSeq[i])...) //translate all the dna sequences to amino acid sequences
	}

	for k := range aaSeq { //for the amino acid(s) that make up each bed record
		for l := range aaSeq[k] { //look through the amino acid sequence for that bed record
			outRecord = append(outRecord, dna.AminoAcidToString(aaSeq[k][l])) //write the three-letter code for that amino acid into the variable to store inal output
		}
	}

	fileio.Write(output, outRecord) //write the contents of the outRecord variable to the file specified in the arguments.

}

// example of usage functions for commands in the gonomics/cmd directory
func usage() {
	fmt.Print(
		"newTool - takes a bed and fasta and converts the bed sequences into amino acide sequences\n" +
			"Usage:\n" +
			"newTool bedFile fastaFile\n" +
			"options:\n")
	flag.PrintDefaults()
}

//example main function for gonomics/cmd directory
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

	bedToAminoAcid(bedFile, fastaFile, outFile)
}
