package main

import (
	"fmt"
	"github.com/vertgenlab/gonomics/fasta"
)

func ReconAccuracy(simFilename string, recFilename string) float64 {
	tot := 0.0
	sim := fasta.Read(simFilename)
	rec := fasta.Read(recFilename)
	des := "descendents_" + simFilename //refers to file created by simulate.RemoveAncestors which has only leaf nodes labelled
	simLeaves := fasta.Read(des)
	for i := 0; i < len(sim); i++ {
		for j := 0; j < len(simLeaves); j++ {
			if sim[i].Name == simLeaves[j].Name { //if current fasta is a leaf fasta
				sim = append(sim[:i], sim[i+1:]...)
				//sim is fastas of every node but the leaf nodes, so we are adding leaf fastas to the original simulation file
			}
		}
	}
	for i := 0; i < len(sim); i++ {
		num := 0.0

		for k := 0; k < len(sim[0].Seq); k++ {
			if sim[i].Seq[k] != rec[i].Seq[k] {
				num = num + 1
			}
		}
		tot = tot + num
	}
	accuracy := tot / (float64(len(sim)) * float64(len(sim[0].Seq))) * 100.0
	acc := 100 - accuracy
	fmt.Print("accuracy over all nodes= ", acc, "%", "\n")
	return acc
}
