package reconstruct

import (
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/expandedTree"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/simulate"
	"testing"
)

var ReconstructTests = []struct {
	NewickFileName       string
	GenePredFile         string
	RandFa               string
	RandFaSeqName        string
	SimTree              string
	LeavesFile           string
	ReconOutFile         string
	GcContent            float64
	Length               int
	BiasLeafName         string
	NonBiasProbThreshold float64
	HighestProbThreshold float64
}{
	{NewickFileName: "testdata/newickLongBranches.txt",
		GenePredFile:         "testdata/genePred.gp",
		RandFa:               "testdata/RandGeneOutput.fa",
		RandFaSeqName:        "test",
		SimTree:              "testdata/simOut.fa",
		LeavesFile:           "testdata/leavesOnly.fa",
		ReconOutFile:         "testdata/reconOut.fa",
		GcContent:            0.42,
		Length:               1005,
		BiasLeafName:         "",
		NonBiasProbThreshold: 0,
		HighestProbThreshold: 0},
}

func TestReconstruct(t *testing.T) {
	var leaves []*expandedTree.ETree
	var tree *expandedTree.ETree
	var err error
	var accuracyData map[string]float64
	var baseAccuracy map[string][]float64
	var baseAccData []float64
	var foundInMap bool

	for _, v := range ReconstructTests {
		tree, err = expandedTree.ReadNewick(v.NewickFileName)
		exception.PanicOnErr(err)
		fasta.Write(v.RandFa, simulate.RandGene(v.RandFaSeqName, v.Length, v.GcContent))
		simulate.Simulate(v.RandFa, tree, v.GenePredFile, false)
		WriteTreeToFasta(tree, v.SimTree)
		WriteLeavesToFasta(tree, v.LeavesFile)

		tree, err = expandedTree.ReadTree(v.NewickFileName, v.LeavesFile)
		exception.FatalOnErr(err)
		leaves = expandedTree.GetLeaves(tree)
		for i := 0; i < len(leaves[0].Fasta.Seq); i++ {
			LoopNodes(tree, i, v.BiasLeafName, v.NonBiasProbThreshold, v.HighestProbThreshold)
		}
		WriteTreeToFasta(tree, v.ReconOutFile)

		accuracyData, baseAccuracy = ReconAccuracy(v.SimTree, v.ReconOutFile, v.LeavesFile, v.GenePredFile, true)
		for name, acc := range accuracyData {
			if name == "D(leaf)" || name == "E(leaf)" || name == "B(leaf)" {
				if acc != 100 {
					t.Errorf("Accuracy for D, E and B should be 100, but accuracy for %s is: %f.", name, acc)
				}
			}
		}

		baseAccData, foundInMap = baseAccuracy["A"]
		if !foundInMap {
			t.Error("node A not found in baseAccuracy data, check tree input.")
		} else if baseAccData[0] < 97.3 || baseAccData[0] > 98.0 {
			t.Errorf("First base accuracy for A in tree should be 97.313433, but is %f.", baseAccData[0])
		}

		baseAccData, foundInMap = baseAccuracy["D"]
		if !foundInMap {
			t.Error("Node D not found in baseAccuracy data, check tree input.")
		} else if baseAccData[0] != 100 {
			t.Errorf("First base accuracy for D should be 100.0, but if %f.", baseAccData[0])
		}

		fileio.EasyRemove(v.RandFa)
		fileio.EasyRemove(v.LeavesFile)
		fileio.EasyRemove(v.ReconOutFile)
		fileio.EasyRemove(v.SimTree)
	}
}
