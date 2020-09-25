package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/reconstruct"
	"log"
)

//TODO: this code will need to change drastically for sequences of varying lengths.
//The loop through the sequence is restricted by a single fasta and the tot calculation will need to calculate the total number of bps
//ReconAccuracy calculates the total number of incorrectly reconstructed base pairs in a tree and returns a percentage of correct base calls
func ReconAccuracy(simFilename string, reconFilename string, outFilename string) {
	answer := reconstruct.ReconAccuracy(simFilename, reconFilename)
	out := fileio.EasyCreate(outFilename)
	defer out.Close()

	for name, accuracy := range answer {
		fmt.Fprintf(out, "%s %f \n", name, accuracy)
	}
}

func usage() {
	fmt.Print(
		"ReconAccuracy takes in a fasta file of simulated evolution along a tree, and the reconstructed fastas of the same tree and returns the percentage accuracy of the sequences of each node and all nodes in the tree.\n" +
			"reconAccuracy <simulationOut.fasta> <reconstructionOut.fasta> <outFilename.txt> \n" +
			"options:\n")
	flag.PrintDefaults()
}

func main() {
	var expectedNumArgs = 3

	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime | log.Lshortfile)
	flag.Parse()

	if len(flag.Args()) != expectedNumArgs {
		flag.Usage()
		log.Fatalf("Error: expecting %d arguments, but got %d\n",
			expectedNumArgs, len(flag.Args()))
	}

	simulationOut := flag.Arg(0)
	reconstructionOut := flag.Arg(1)
	outFile := flag.Arg(2)

	ReconAccuracy(simulationOut, reconstructionOut, outFile)
}
