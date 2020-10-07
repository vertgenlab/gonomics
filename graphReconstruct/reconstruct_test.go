package graphReconstruct

import (
	"github.com/vertgenlab/gonomics/expandedTree"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/simulate"
	"log"
	"testing"
)

var GCcontent = 0.42
var input = []struct {
	newickFilename string // first input
	length         int    // second input
}{
	{"newTestFiles/newickLongBranches.txt", 1005},
}

func Test_reconstruct(t *testing.T) {
	for _, test := range input {
		tre, er := expandedTree.ReadNewick(test.newickFilename)
		if er != nil {
			log.Fatal("Couldn't read file")
		}
		fasta.Write("RandGeneOutput.fasta", simulate.RandGene("test", test.length, GCcontent)) //galGal6 GC
		simulate.Simulate("RandGeneOutput.fasta", tre)
		WriteTreeToFasta(tre, "simOut.fasta")
		WriteLeavesToFasta(tre, "leavesOnly.Fasta")

		tr := expandedTree.ReadTree(test.newickFilename, "leavesOnly.fasta")
		leaves := expandedTree.GetLeaves(tr)
		for i := 0; i < len(leaves[0].Fasta.Seq); i++ {
			LoopNodes(tr, i)
		}
		WriteTreeToFasta(tr, "reconOut.fasta")
		accuracyData := ReconAccuracy("simOut.fasta", "reconOut.fasta")
		for name, accuracy := range accuracyData {
			log.Printf("%s %f \n", name, accuracy)
		}
	}
}
