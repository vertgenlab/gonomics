package reconstructSeq

import (
	"github.com/vertgenlab/gonomics/reconstruct"
	"github.com/vertgenlab/gonomics/expandedTree"
	"github.com/vertgenlab/gonomics/fasta"
)

func ReconstructSeq (newickInput string, fastaInput string, outputFilename string) {
	tree := expandedTree.ReadTree(newickInput, fastaInput)
	leaves := expandedTree.GetLeaves(tree)
	branches := expandedTree.GetBranch(tree)
	var treeFastas []*fasta.Fasta
	
	for i := 0; i < len(leaves[0].Fasta.Seq); i++ {
		reconstruct.LoopNodes(tree, i)
	}
	for j := 0; j < len(leaves); j++ {
		treeFastas = append(treeFastas, leaves[j].Fasta)
	}
	for k := 0; k < len(branches); k++ {
		treeFastas = append(treeFastas, branches[k].Fasta)
	}
	fasta.Write(outputFilename, treeFastas)
}

func Accuracy (simFilename string, reconFilename string) 