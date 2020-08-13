package reconstruct

import (
	"github.com/vertgenlab/gonomics/simulate"
	"github.com/vertgenlab/gonomics/tree_newick"
	"testing"
)

var input = []struct {
	newick_filename string // first input
	length          int    // second input
}{
	{"testdata/test_newick", 33},
	{"testdata/test_newick", 333},
	{"testdata/hackett_newick", 66}}

func Test_reconstruct(t *testing.T) {
	for _, test := range input {

		//example of how to run simulation: 1) read in tree (no fastas) 2) simulate random gene 3) simulate evolution 4) remove ancestors for reconstruction 4) set up fastas
		tre, er := tree_newick.ReadNewick(test.newick_filename)
		if er != nil {
		}
		simulate.Rand_gene("test", test.length)
		simulate.Simulate("test.fasta", "test_tree.fasta", tre)
		simulate.Remove_ancestors("test_tree.fasta", tre)
		tree_newick.Set_fastas_up(tre, "test_tree.fasta")

		//example of how to run reconstruction: 1) reconstruct 2) check accuracy
		tr := tree_newick.Read_tree(test.newick_filename, "test_tree.fasta")
		Reconstruct(tr, "test_reconstruction.fasta")

		//example of check simulation vs reconstrucion
		a := Accuracy("test_tree.fasta", "test_reconstruction.fasta")

		if a < 0.8 {
			t.Errorf("Accuracy of %v is below 0.8, expected is (~ .9, .91, .97)", a)
		}

	}
}
