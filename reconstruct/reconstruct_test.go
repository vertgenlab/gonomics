package reconstruct

import (
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/expandedTree"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/simulate"
	"log"
	"testing"
)

var GCcontent = 0.42
var input = []struct {
	newickFilename string // first input
	length         int    // second input
}{
	{"testdata/newickLongBranches.txt", 1005},
}

func Test_reconstruct(t *testing.T) {
	var leaves []*expandedTree.ETree
	var accuracyData map[string]float64
	for _, test := range input {
		tre, er := expandedTree.ReadNewick(test.newickFilename)
		if er != nil {
			log.Fatal("Couldn't read file")
		}
		fasta.Write("RandGeneOutput.fasta", simulate.RandGene("test", test.length, GCcontent)) //galGal6 GC
		simulate.Simulate("RandGeneOutput.fasta", tre, "testdata/genePred.gp", false)
		WriteTreeToFasta(tre, "simOut.fasta")
		WriteLeavesToFasta(tre, "leavesOnly.fasta")

		tr, err := expandedTree.ReadTree(test.newickFilename, "leavesOnly.fasta")
		exception.FatalOnErr(err)
		leaves = expandedTree.GetLeaves(tr)
		for i := 0; i < len(leaves[0].Fasta.Seq); i++ {
			LoopNodes(tr, i)
		}
		WriteTreeToFasta(tr, "reconOut.fasta")
		accuracyData, _ = ReconAccuracy("simOut.fasta", "reconOut.fasta", "leavesOnly.fasta", "testdata/genePred.gp", false)
		for name, accuracy := range accuracyData {
			log.Printf("%s %f \n", name, accuracy)
		}
	}

	if accuracyData[leaves[0].Name] != 100 {

	}

	fileio.EasyRemove("RandGeneOutput.fasta")
}

func TestReconAccuracyByBase(t *testing.T) {
	_, baseAccuracy := ReconAccuracy("simOut.fasta", "reconOut.fasta", "leavesOnly.fasta", "testdata/genePred.gp", true)

	for name, data := range baseAccuracy {
		for d := range data {
			if d == 0 {
				log.Printf("%s First Base %f \n", name, data[d])
			} else if d == 1 {
				log.Printf("%s Second Base %f \n", name, data[d])
			} else {
				log.Printf("%s Third Base %f \n", name, data[d])
			}
		}
	}

	fileio.EasyRemove("leavesOnly.fasta")
	fileio.EasyRemove("reconOut.fasta")
	fileio.EasyRemove("simOut.fasta")
}
