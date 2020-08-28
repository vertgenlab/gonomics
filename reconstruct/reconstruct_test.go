package reconstruct

import (
	"github.com/vertgenlab/gonomics/expandedTree"
	"github.com/vertgenlab/gonomics/fasta"
	//"github.com/vertgenlab/gonomics/simulate"
	"log"
	"testing"
)

var GCcontent = 0.42
var input = []struct {
	newick_filename string // first input
	length          int    // second input
}{
	//{"testdata/test_newick", 33},
	{"testdata/test_newick", 333},
	{"testdata/hackett_newick", 66},
}

func TestTinyRecon(t *testing.T) {
	//example of how to run simulation: 1) read in tree (no fastas) 2) simulate random gene 3) simulate evolution 4) remove ancestors for reconstruction 4) assign fastas
	log.Printf("Reading in files...\n")
	tr := expandedTree.ReadTree("testdata/tinyApe.nh", "testdata/tinyApe.fa")

	leaves := expandedTree.GetLeaves(tr)
	branches := expandedTree.GetBranch(tr)
	var treeFastas []*fasta.Fasta

	seqLength := len(leaves[0].Fasta.Seq)
	for i := 0; i < seqLength; i++ {
		log.Printf("Setting leaf state for column %d...\n", i)
		SetLeafState(tr, i)
		log.Printf("Setting internal state for column %d...\n", i)
		SetInternalState(tr)
		log.Printf("Looping over nodes for column %d...\n", i)
		LoopNodes(tr)
	}
	for j := 0; j < len(leaves); j++ {
		treeFastas = append(treeFastas, leaves[j].Fasta)
	}
	for j := 0; j < len(branches); j++ {
		treeFastas = append(treeFastas, branches[j].Fasta)
	}
	fasta.Write("testdata/tinyApeInternal.fa", treeFastas)
}

/* commenting this out to do an easy test to start
func Test_reconstruct(t *testing.T) {
	for _, test := range input {

		//example of how to run simulation: 1) read in tree (no fastas) 2) simulate random gene 3) simulate evolution 4) remove ancestors for reconstruction 4) assign fastas
		tre, er := expandedTree.ReadNewick(test.newick_filename)
		common.ExitIfError(er)
		fasta.Write("test.fasta", simulate.RandGene("test", test.length, GCcontent)) //galGal6 GC
		simulate.Simulate("test.fasta", "test_tree.fasta", tre)
		simulate.RemoveAncestors("test_tree.fasta", tre)
		expandedTree.AssignFastas(tre, "test_tree.fasta")

		//example of how to run reconstruction: 1) reconstruct 2) check accuracy
		tr := expandedTree.ReadTree(test.newick_filename, "test_tree.fasta")
		LoopNodes(tr)

		//example of check simulation vs reconstrucion
		a := Accuracy("test_tree.fasta", "test_reconstruction.fasta")

		if a < 0.8 {
			t.Errorf("Accuracy of %v is below 0.8, expected is (~ .9, .91, .97)", a)
		}

	}
}*/
