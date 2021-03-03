package expandedTree

import (
	"log"
	"testing"
)

var actualName = "root"

func TestReadTree(t *testing.T) {
	tree, err := ReadTree("testdata/newick.txt", "testdata/fasta.fasta")
	if err != nil {
		log.Print("Problem in ReadTree")
	} else {
		if tree.Name != actualName {
			log.Fatal("Incorrect name in tree")
		}
		if tree.Name == "A" {
			if tree.BranchLength != 0.10 {
				log.Fatal("Wrong branch length for leaf A")
			}
		}
	}
}
