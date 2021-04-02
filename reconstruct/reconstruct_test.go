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
	var name string
	var data []float64

	for name, data = range baseAccuracy {
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

	data, ok := baseAccuracy["A"]
	if !ok {
		t.Error("node A not found in baseAccuracy data, check tree input.")
	} else {
		if data[0] < 97.0 {
			t.Errorf("First base accuracy for A in tree should be 97.910448, but is %f.", data[0])
		}
	}
	data, ok = baseAccuracy["D"]
	if !ok {
		t.Error("Node D not found in baseAccuracy data, check tree input.")
	} else {
		if data[0] != 100 {
			t.Errorf("First base accuracy for D should be 100.0, but if %f.", data[0])
		}
	}

	fileio.EasyRemove("leavesOnly.fasta")
	fileio.EasyRemove("reconOut.fasta")
	fileio.EasyRemove("simOut.fasta")
}
