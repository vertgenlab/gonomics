// Converts a bed and fasta sequences into amino acid sequences within the bed regions

// Command Group: "Data Conversion"

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

// bedToAminoAcid is a complete function made from gonomics packages, it takes a bed and a fasta and returns the amino acid sequence of all the beds together.
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
			outRecord = append(outRecord, dna.AminoAcidToString(aaSeq[k][l])) //write the three-letter code for that amino acid into the variable to store in output
		}
	}

	fileio.Write(output, outRecord) //write the contents of the outRecord variable to the file specified in the arguments.

}

// example of usage functions for commands in the gonomics/cmd directory
func usage() {
	fmt.Print(
		"bedToAminoAcid - takes a bed and fasta and converts the bed sequences into amino acid sequences\n" +
			"Usage:\n" +
			"bedToAminoAcid bedFile fastaFile\n" +
			"options:\n")
	flag.PrintDefaults()
}

// example main function for gonomics/cmd directory
func main() {
	var expectedNumArgs int = 3 //this line specifies the expected number of arguments for the function. Any deviation from this (not including options) will produce an error that prints the usage statement

	flag.Usage = usage                                   //if a usage flag is thrown (as shown in the if statement below) the usage statement will be printed to screen
	log.SetFlags(log.Ldate | log.Ltime | log.Lshortfile) //this specifies that any log package output must include anything set in the arguments of this function
	flag.Parse()                                         //this parses the different arguments given to command

	if len(flag.Args()) != expectedNumArgs { //this is the if statement that catches if the number of arguments is wrong and prints the usage
		flag.Usage() //print usage statement
		log.Fatalf("Error: expecting %d arguments, but got %d\n",
			expectedNumArgs, len(flag.Args())) //fail and print statement
	}

	bedFile := flag.Arg(0) // the numbered ordered of these variables is to tell the program how to order the inputs, the first will be used in the function as the bedFile for example
	fastaFile := flag.Arg(1)
	outFile := flag.Arg(2)

	bedToAminoAcid(bedFile, fastaFile, outFile) //this is where the function is called so that the program runs if there is now error beforehand. Without this, nothing would run and the function would be flagged as being declared but not used
}
