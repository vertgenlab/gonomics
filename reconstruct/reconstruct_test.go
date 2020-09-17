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
		simulate.Simulate("RandGeneOutput.fasta", tre)                                         //produced whole tree fasta file
		simulate.RemoveAncestors("simOut.fasta", tre, "leavesOnly.fasta")                      //removes everything but leaves from simOut and produced another file
		expandedTree.AssignFastas(tre, "simOut.fasta")
		//TODO: none of this works this way anymore, should call simulateEvol command to run these.

		//example of how to run reconstruction: 1) reconstruct 2) check accuracy
		tr := expandedTree.ReadTree(test.newickFilename, "simOut.fasta")
		for i := 0; i < len(tr.Fasta.Seq); i++ {
			LoopNodes(tr, i)
		}
	}
}
