package main

import (
	"fmt"
	"github.com/vertgenlab/gonomics/fasta"
)

func ReconAccuracy(simFilename string, reconFilename string) float64 {
	tot := 0.0
	sim := fasta.Read(simFilename)
	recon := fasta.Read(reconFilename)
	//TODO: add usages
	for i := 0; i < len(sim); i++ {
		num := 0.0

		for k := 0; k < len(sim[0].Seq); k++ {
			if sim[i].Seq[k] != recon[i].Seq[k] {
				num = num + 1
			}
		}
		tot = tot + num
	}
	accuracy := tot / (float64(len(sim)) * float64(len(sim[0].Seq))) * 100.0
	acc := 100 - accuracy
	fmt.Print("accuracy over all nodes= ", acc, "%", "\n")
	return acc //probably a better way to do it than to print a number to the screen?
}
