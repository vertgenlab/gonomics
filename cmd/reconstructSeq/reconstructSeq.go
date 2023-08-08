// Command Group: "Sequence Evolution & Reconstruction"

// Reconstruct ancient sequences using extant genomes and a newick tree with branch lengths
package main

import (
	"flag"
	"fmt"
	"log"

	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/expandedTree"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/reconstruct"
)

type Settings struct {
	NewickInput          string
	FastaInput           string
	OutFile              string
	BiasLeafName         string
	NonBiasProbThreshold float64
	HighestProbThreshold float64
}

func ReconstructSeq(s Settings) {
	var treeFastas []fasta.Fasta

	if s.NonBiasProbThreshold < 0 || s.NonBiasProbThreshold > 1 {
		log.Fatalf("Error: nonBiasProbThreshold must be a value between 0 and 1. Found: %v.\n", s.NonBiasProbThreshold)
	}

	if s.NonBiasProbThreshold > 0 && s.BiasLeafName == "" {
		log.Fatalf("Error: nonBiasProbThreshold was set, but no BiasLeafName was provided.\n")
	}

	if s.HighestProbThreshold < 0 || s.HighestProbThreshold > 1 {
		log.Fatalf("Error: highestProbThreshold must be a value between 0 and 1. Found: %v.\n", s.HighestProbThreshold)
	}

	tree, err := expandedTree.ReadTree(s.NewickInput, s.FastaInput)
	exception.FatalOnErr(err)

	leaves := expandedTree.GetLeaves(tree)
	branches := expandedTree.GetBranch(tree)

	for i := range leaves[0].Fasta.Seq {
		reconstruct.LoopNodes(tree, i, s.BiasLeafName, s.NonBiasProbThreshold, s.HighestProbThreshold)
	}
	for j := range leaves {
		treeFastas = append(treeFastas, *leaves[j].Fasta)
	}
	for k := range branches {
		treeFastas = append(treeFastas, *branches[k].Fasta)
	}
	fasta.Write(s.OutFile, treeFastas)
}

func usage() {
	fmt.Print(
		"reconstructSeq performs ancestral sequence reconstruction based on an input multiFa alignment and Newick tree." +
			"This program returns a fasta file containing sequences for all nodes of the tree, including the input sequences (leaves)," +
			"and the inferred ancestral nodes.\n" +
			"Usage:\n" +
			"reconstructSeq <treeFile.txt> <inFasta.fasta> <outFasta.fasta>\n" +
			"options:\n")
	flag.PrintDefaults()
}

func main() {
	var biasLeafName *string = flag.String("biasLeafName", "", "Specify an extant (leaf) sequence towards which we will bias reconstruction for the immediate ancestor (parent node).")
	var nonBiasProbThreshold *float64 = flag.Float64("nonBiasBaseThreshold", 0, "When biasing reconstruction for a leaf, this parameter sets the sum of probabilities for non-leaf bases. TODO: make this more clear")
	var highestProbThreshold *float64 = flag.Float64("highestProbThreshold", 0, "The highest probability base must be above this value to be accepted for reconstruction. Otherwise, dna.N will be returned.")

	var expectedNumArgs = 3

	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime | log.Lshortfile)
	flag.Parse()

	if len(flag.Args()) != expectedNumArgs {
		flag.Usage()
		log.Fatalf("Error: expecting %d arguments, but got %d\n",
			expectedNumArgs, len(flag.Args()))
	}

	newickInput := flag.Arg(0)
	fastaInput := flag.Arg(1)
	outFile := flag.Arg(2)

	s := Settings{
		NewickInput:          newickInput,
		FastaInput:           fastaInput,
		OutFile:              outFile,
		BiasLeafName:         *biasLeafName,
		NonBiasProbThreshold: *nonBiasProbThreshold,
		HighestProbThreshold: *highestProbThreshold,
	}

	ReconstructSeq(s)
}
