package reconstruct

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
	//{"newTestFiles/newickShortBranches.txt", 333},
	//{"newTestFiles/newickShortBranches.txt", 66}
}

func Test_reconstruct(t *testing.T) {
	for _, test := range input {

		//example of how to run simulation: 1) read in tree (no fastas) 2) simulate random gene 3) simulate evolution 4) remove ancestors for reconstruction 4) assign fastas
		tre, er := expandedTree.ReadNewick(test.newickFilename)
		if er != nil {
			log.Fatal("Couldn't read file")
		}
		fasta.Write("RandGeneOutput.fasta", simulate.RandGene("test", test.length, GCcontent)) //galGal6 GC
		simulate.Simulate("RandGeneOutput.fasta", tre) //produced whole tree fasta file
		WriteTreeToFasta(tre, "simOut.fasta")
		WriteLeavesToFasta(tre, "leavesOnly.Fasta")

		tr := expandedTree.ReadTree(test.newickFilename, "leavesOnly.fasta")
		leaves := expandedTree.GetLeaves(tr)
		for i := 0; i < len(leaves[0].Fasta.Seq); i++ {
			LoopNodes(tr, i)
		}
		WriteTreeToFasta(tr, "reconOut.fasta")
		log.Printf("accuracy for all nodes: %f percent", ReconAccuracy("simOut.fasta", "reconOut.fasta"))
		accuracyData := ReconAccuracy("simOut.fasta", "reconOut.fasta")
		for name, accuracy := range accuracyData {
		log.Printf("%s %f \n", name, accuracy)
		}
	}
}