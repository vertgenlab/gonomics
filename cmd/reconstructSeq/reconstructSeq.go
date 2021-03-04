package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/expandedTree"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/reconstruct"
	"log"
)

func ReconstructSeq(newickInput string, fastaInput string, outputFilename string) {
	tree, err := expandedTree.ReadTree(newickInput, fastaInput)
	if err != nil {
		log.Fatalf("Error in ReadTree: %e", err)
	}
	leaves := expandedTree.GetLeaves(tree)
	branches := expandedTree.GetBranch(tree)
	var treeFastas []fasta.Fasta

	for i := 0; i < len(leaves[0].Fasta.Seq); i++ {
		reconstruct.LoopNodes(tree, i)
	}
	for j := 0; j < len(leaves); j++ {
		treeFastas = append(treeFastas, *leaves[j].Fasta)
	}
	for k := 0; k < len(branches); k++ {
		treeFastas = append(treeFastas, *branches[k].Fasta)
	}
	fasta.Write(outputFilename, treeFastas)
}

func usage() {
	fmt.Print(
		"reconstructSeq takes in a newick tree file and a fasta file containing fastas for all leaf nodes and returns a fasta file containing fastas for all nodes of the tree\n" +
			"Usage:\n" +
			"reconstructSeq <treeFile.txt> <inFasta.fasta> <outFasta.fasta>\n" +
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

	treeFile := flag.Arg(0)
	fastaFile := flag.Arg(1)
	outFile := flag.Arg(2)

	ReconstructSeq(treeFile, fastaFile, outFile)
}
