package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/expandedTree"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/simulate"
	"log"
)

//TODO: option for seeded or unseeded random numbers (include in simulate.go)
func SimulateEvol(rootFastaFile string, treeFile string, gp string, simOutFile string, leafOutFile string) {
	tree, err := expandedTree.ReadTree(treeFile, rootFastaFile)
	if err != nil {
		log.Fatalf("Error in ReadTree: %e", err)
	}
	var fastas []fasta.Fasta
	var leafFastas []fasta.Fasta
	simulate.Simulate(rootFastaFile, tree, gp, true)
	nodes := expandedTree.GetTree(tree)

	for i := 0; i < len(nodes); i++ {
		fastas = append(fastas, *nodes[i].Fasta)
		if nodes[i].Left == nil && nodes[i].Right == nil {
			leafFastas = append(leafFastas, *nodes[i].Fasta)
		}
	}
	fasta.Write(simOutFile, fastas)
	fasta.Write(leafOutFile, leafFastas)
}

func usage() {
	fmt.Print(
		"simulateEvol takes in a root fasta, the sequence's genePred file, and a newick formatted tree with branch lengths and simulates evolution along the tree. It returns a list of fastas for the whole tree for reference and a list of fastas from leaves for reconstruction.\n" +
			"Usage:\n" +
			"simulateEvol <rootFasta.fasta> <newickTree.txt> <genePred filename> <outFile.fasta> <leafOutputFile.fasta> \n" +
			"options:\n")
	flag.PrintDefaults()
}

func main() {
	var expectedNumArgs = 5

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
	gene := flag.Arg(2)
	outFile := flag.Arg(3)
	leafOutFile := flag.Arg(4)

	SimulateEvol(rootFasta, newickTree, gene, outFile, leafOutFile)
}
