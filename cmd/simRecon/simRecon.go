package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/expandedTree"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/reconstruct"
	"github.com/vertgenlab/gonomics/simulate"
	"log"
)

//SimulateEvolve takes in a root fasta file, a newick tree, and a genePred file and writes two files, one with simulated sequences for every node of the tree
//and one for simulated sequences on just the leaves of the tree to be used for reconstruction
func SimulateEvolve(rootFastaFile string, treeFile string, gp string, simOutFile string, leafOutFile string) {
	tree, err := expandedTree.ReadTree(treeFile, rootFastaFile)
	exception.FatalOnErr(err)
	var fastas []fasta.Fasta
	var leafFastas []fasta.Fasta
	simulate.Simulate(rootFastaFile, tree, gp, false)
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

//ReconstructSeq takes in a newick tree and leaf sequences (output of simulateEvolve) and writes a file containing reconstructed sequences for the whole input tree
func ReconstructSeq(newickInput string, fastaInput string, outputFilename string) {
	tree, err := expandedTree.ReadTree(newickInput, fastaInput)
	exception.FatalOnErr(err)
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

//SimRecon simulates evolution, performs reconstruction, and then evaluates the accuracy of the reconstruction
//output for accuracy is a file with 3 accuracy calculations for each ancestor: exon-specific, non-exon-specific and a total for that ancestor
//this code is designed for sequences which cannot change in length as they are simulated or reconstructed
func SimRecon(rootFastaFile string, treeFile string, gp string, simOutFile string, leafOutFile string, reconOutFile string, accuracyOutFile string) {
	SimulateEvolve(rootFastaFile, treeFile, gp, simOutFile, leafOutFile)
	ReconstructSeq(treeFile, leafOutFile, reconOutFile)

	answer := reconstruct.ReconAccuracy(simOutFile, reconOutFile, leafOutFile, gp)
	out := fileio.EasyCreate(accuracyOutFile)

	for name, accuracy := range answer {
		fmt.Fprintf(out, "%s\t%f\n", name, accuracy)
	}
	err := out.Close()
	exception.PanicOnErr(err)
}

func usage() {
	fmt.Print(
		"SimRecon takes in a root fasta file and a tree. It simulates evolution along the given tree, performs reconstruction based on the leaf nodes of the tree and then calculates accuracy at each node. All outputs are written to a file.\n" +
			"simRecon <rootFile.fasta> <treeFile.txt> <genePred.gp> <simOutFile.fasta> <leafOutFile.fasta> <reconOutFile.fasta> <accuracyOutFile.txt> \n" +
			"options:\n")
	flag.PrintDefaults()
}

func main() {
	var expectedNumArgs = 7

	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime | log.Lshortfile)
	flag.Parse()

	if len(flag.Args()) != expectedNumArgs {
		flag.Usage()
		log.Fatalf("Error: expecting %d arguments, but got %d\n",
			expectedNumArgs, len(flag.Args()))
	}

	rootFastaFile := flag.Arg(0)
	treeFile := flag.Arg(1)
	gp := flag.Arg(2)
	simOutFile := flag.Arg(3)
	leafOutFile := flag.Arg(4)
	reconOutFile := flag.Arg(5)
	accuracyOutFile := flag.Arg(6)

	SimRecon(rootFastaFile, treeFile, gp, simOutFile, leafOutFile, reconOutFile, accuracyOutFile)
}
