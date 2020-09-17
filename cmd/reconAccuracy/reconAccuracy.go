package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/fasta"
	"log"
)

func ReconAccuracy(simFilename string, reconFilename string) float64 {
	//TODO: print to file?
	tot := 0.0
	sim := fasta.Read(simFilename)
	recon := fasta.Read(reconFilename)
	//TODO: account for inflating accuracy by counting leaf nodes?
	//TODO: leaf nodes are not the same between sim and recon files
	for i := 0; i < len(sim); i++ {
		num := 0.0
		for k := 0; k < len(sim[i].Seq); k++ { //should that be sim[0] or sim[i]
			if sim[i].Name == recon[i].Name { //added this line bc i noticed the output from recon was in a different order than output from simulation
				if sim[i].Seq[k] != recon[i].Seq[k] {
					num = num + 1
				}
			}
		}
		tot = tot + num
	}
	accuracy := tot / (float64(len(sim)) * float64(len(sim[0].Seq))) * 100.0
	fmt.Printf("tot: %f, len(sim): %f, len(sim[0].Seq): %f \n", tot, float64(len(sim)), float64(len(sim[0].Seq)))
	acc := 100 - accuracy
	fmt.Print("accuracy over all nodes= ", acc, "%", "\n")
	return acc //probably a better way to do it than to print a number to the screen?
}

func usage() {
	fmt.Print(
		"ReconAccuracy takes in a fasta file of simulated evolution along a tree, and the reconstructed fastas of the same tree and returns the percentage accuracy of the sequences of the tree.\n" +
			"reconAccuracy <simulationOut.fasta> <reconstructionOut.fasta> \n" +
			//TODO: update if printing to another file
			"options:\n")
	flag.PrintDefaults()
}

func main() {
	var expectedNumArgs = 2

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

	ReconAccuracy(simulationOut, reconstructionOut)
}
