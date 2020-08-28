package reconstruct

import (
	"github.com/vertgenlab/gonomics/expandedTree"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/simulate"
	"testing"
)

var GCcontent = 0.42
var input = []struct {
	newick_filename string // first input
	length          int    // second input
}{
	{"testdata/test_newick", 33},
	{"testdata/test_newick", 333},
	{"testdata/hackett_newick", 66}}

func Test_reconstruct(t *testing.T) {
	for _, test := range input {

		//example of how to run simulation: 1) read in tree (no fastas) 2) simulate random gene 3) simulate evolution 4) remove ancestors for reconstruction 4) assign fastas
		tre, er := expandedTree.ReadNewick(test.newick_filename)
		if er != nil {
		}
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
}
