package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/expandedTree"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/simulate"
	"log"
)

func SimulateEvol(rootFastaFile string, treeFile string, simOutFile string, leafOutFile string) {
	tree := expandedTree.ReadTree(treeFile, rootFastaFile)
	var fastas []*fasta.Fasta
	simulate.Simulate(rootFastaFile, tree)
	nodes := expandedTree.GetTree(tree)

	for i := 0; i < len(nodes); i++ {
		fastas = append(fastas, nodes[i].Fasta)
	}
	fasta.Write(simOutFile, fastas)

	simulate.RemoveAncestors(simOutFile, tree, leafOutFile)
}

func usage() {
	fmt.Print(
		"simulateEvol takes in a root fasta and a newick formatted tree with branch lengths and simulates evolution along the tree. It returns a list of fastas for the whole tree for reference and a list of fastas from leaves for reconstruction.\n" +
			"Usage:\n" +
			"simulateEvol <rootFasta.fasta> <newickTree.txt> <outFile.fasta> <leafOutputFile.fasta> \n" +
			"options:\n")
	flag.PrintDefaults()
}

func main() {
	var expectedNumArgs = 4

	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime | log.Lshortfile)
	flag.Parse()

	if len(flag.Args()) != expectedNumArgs {
		flag.Usage()
		log.Fatalf("Error: expecting %d arguments, but got %d\n",
			expectedNumArgs, len(flag.Args()))
	}

	rootFasta := flag.Arg(0)
	newickTree := flag.Arg(1)
	outFile := flag.Arg(2)
	leafOutFile := flag.Arg(3)

	SimulateEvol(rootFasta, newickTree, outFile, leafOutFile)
}
