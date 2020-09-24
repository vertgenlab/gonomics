package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/fasta"
	"log"
)

//TODO: this code will need to change drastically for sequences of varying lengths.
//The loop through the sequence is restricted by a single fasta and the tot calculation will need to calculate the total number of bps
//ReconAccuracy calculates the total number of incorrectly reconstructed base pairs in a tree and returns a percentage of correct base calls
func ReconAccuracy(simFilename string, reconFilename string) {
	var allNodes string
	allNodes = "all Nodes"
	var found bool = false
	tot := 0.0
	sim := fasta.Read(simFilename)
	recon := fasta.Read(reconFilename)

	answer := make(map[string]float64)

	for i := 0; i < len(sim); i++ {
		num := 0.0
		found = false
		for j := 0; j < len(recon); j++ {
			if sim[i].Name == recon[j].Name {
				found = true
				//DEBUG: log.Printf("\n%s \n%s \n", dna.BasesToString(sim[i].Seq), dna.BasesToString(recon[j].Seq))
				for k := 0; k < len(sim[0].Seq); k++ {
					if sim[i].Seq[k] != recon[j].Seq[k] {
						num = num + 1
					}
				}
			}
		}
		if found == false {
			log.Fatal("Did not find all simulated sequences in reconstructed fasta.")
		} else {
			accuracy := num / float64(len(sim[i].Seq)) * 100.0
			//DEBUG: fmt.Printf("tot: %f, len(sim): %f, len(sim[0].Seq): %f \n", tot, float64(len(sim)), float64(len(sim[0].Seq)))
			acc := 100 - accuracy
			answer[sim[i].Name] = acc
		}
		tot = tot + num
	}
	accuracy := tot / (float64(len(sim)) * float64(len(sim[0].Seq))) * 100.0
	//DEBUG: fmt.Printf("tot: %f, len(sim): %f, len(sim[0].Seq): %f \n", tot, float64(len(sim)), float64(len(sim[0].Seq)))
	acc := 100 - accuracy
	answer[allNodes] = acc
	for name, accuracy := range answer {
		log.Printf("%s %f \n", name, accuracy)
	}
}

func usage() {
	fmt.Print(
		"ReconAccuracy takes in a fasta file of simulated evolution along a tree, and the reconstructed fastas of the same tree and returns the percentage accuracy of the sequences of each node and all nodes in the tree.\n" +
			"reconAccuracy <simulationOut.fasta> <reconstructionOut.fasta> \n" +
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
