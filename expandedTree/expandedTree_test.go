package expandedTree

import (
	"log"
	"testing"
)

var actualName = "root"

func TestReadTree(t *testing.T) {
	tree, err := ReadTree("testdata/newick.txt", "testdata/fasta.fasta")
	treeNodes := GetTree(tree)
	if err != nil {
		log.Fatal("Problem in ReadTree")
	} else {
		if tree.Name != actualName {
			log.Fatal("Incorrect name in tree")
		}
		for b := 0; b < len(treeNodes); b++ {
			if treeNodes[b].Name == "A" {
				if treeNodes[b].BranchLength != 0.20 {
					log.Fatal("incorrect branch length assignment")
				}
			}
		}
	}
}
