package expandedTree

import (
	"github.com/vertgenlab/gonomics/exception"
	"log"
	"testing"
)

func TestReadTree(t *testing.T) {
	var actualName = "root"
	var aFound = false
	tree, err := ReadTree("testdata/newick.txt", "testdata/fasta.fasta")
	treeNodes := GetTree(tree)
	exception.FatalOnErr(err)
	if tree.Name != actualName {
		t.Error("Incorrect name at root")
	}
	for b := 0; b < len(treeNodes); b++ {
		if treeNodes[b].Name == "A" {
			aFound = true
			if treeNodes[b].BranchLength != 0.20 {
				t.Error("incorrect branch length assignment")
			}
		}
	}
	if !aFound {
		log.Fatal("A node not found in tree")
	}
}
