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

//SimulateEvolve takes in a root fasta file, a newick tree, and gene structure genePred file for the fasta and returns a full simulated tree and a tree with sequence only at the leaves for reconstruction
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

//ReconstructSeq takes in a newick tree and leaf sequences and returns a reconstructed tree
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
func SimRecon(rootFastaFile string, treeFile string, gp string, simOutFile string, leafOutFile string, reconOutFile string, accuracyOutFile string, option int, baseAccFile string) {
	//TODO: this code will need to change drastically for sequences of varying lengths.
	//The loop through the sequence is restricted by the length of a single fasta and the tot calculation will need to calculate the total number of bps
	//ReconAccuracy calculates the total number of incorrectly reconstructed base pairs in a tree and returns a percentage of correct base calls
	SimulateEvolve(rootFastaFile, treeFile, gp, simOutFile, leafOutFile)
	ReconstructSeq(treeFile, leafOutFile, reconOutFile)

	answer, byBaseAnswer := reconstruct.ReconAccuracy(simOutFile, reconOutFile, leafOutFile, gp, option)
	out := fileio.EasyCreate(accuracyOutFile)
	defer out.Close()

	for name, accuracy := range answer {
		fmt.Fprintf(out, "%s\t%f\n", name, accuracy)
	}
	if option == 1 && baseAccFile != "" {
		baseAccOut := fileio.EasyCreate(baseAccFile)
		defer baseAccOut.Close()
		for species, baseAcc := range byBaseAnswer {
			firstBase := species + " First Base"
			secondBase := species + " Second Base"
			thirdBase := species + " Third Base"
			for b, _ := range baseAcc {
				if b == 0 {
					fmt.Fprintf(baseAccOut, "%s\t%f\n", firstBase, baseAcc[b])
				} else if b == 1 {
					fmt.Fprintf(baseAccOut, "%s\t%f\n", secondBase, baseAcc[b])
				} else {
					fmt.Fprintf(baseAccOut, "%s\t%f\n", thirdBase, baseAcc[b])
				}
			}
		}
	}
}

func usage() {
	fmt.Print(
		"SimRecon takes in a root fasta file, a genePred, and a tree. It simulates evolution along the tree, performs reconstruction based on the leaf nodes of the tree and then calculates accuracy at each node.\n" +
			"simRecon [-option] <rootFile.fasta> <treeFile.txt> <genePred.gp> <simOutFile.fasta> <leafOutFile.fasta> <reconOutFile.fasta> <accuracyOutFile.txt>\n" +
			"options:\n")
	flag.PrintDefaults()
}

func main() {
	var expectedNumArgs = 7
	var option = flag.Int("option", 0, "Specify to 1 to return an accuracy breakdown by base position in codons. Default value is 0, which returns normal output.")
	var baseAccFile = flag.String("baseAccFile", "", "Specify a filename for the output of accuracy by position of base in codon (if option = 1).")

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

	SimRecon(rootFastaFile, treeFile, gp, simOutFile, leafOutFile, reconOutFile, accuracyOutFile, *option, *baseAccFile)
}
