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
	SubMatrix            bool
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
		HighestProbThreshold: 0,
		SubMatrix:            false,
	},
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
			LoopNodes(tree, i, v.BiasLeafName, v.NonBiasProbThreshold, v.HighestProbThreshold, v.SubMatrix)
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

// this benchmark takes about 60s CPU to run, so I've commented it out. Go ahead
// and uncomment if you're interested in reproducing it or altering it. -RJM
/*
var EmpiricalReconstructionComparison = []struct {
	TestName                    string
	SimSubstitutionMatrixFile   string
	ReconSubstitutionMatrixFile string
	LeavesFile                  string
	NewickFile                  string
	ReconOutFile                string
	NumTrees                    int
	NodeGammaAlpha              float64
	NodeGammaBeta               float64
	BranchAlpha                 float64
	BranchBeta                  float64
	SetSeed                     int64
	GcContent                   float64
	BiasLeafName                string
	NonBiasProbThreshold        float64
	HighestProbThreshold        float64
	SeqLen                      int
	SubMatrix                   bool
}{
	{TestName: "GtrSimJukesRecon",
		SimSubstitutionMatrixFile:   "testdata/substitutionMatrix.txt",
		ReconSubstitutionMatrixFile: "testdata/jukesCantor.txt",
		LeavesFile:                  "testdata/leavesFile.txt",
		NewickFile:                  "testdata/currNewick.txt",
		ReconOutFile:                "testdata/reconOutFile.txt",
		NumTrees:                    50,
		NodeGammaAlpha:              3,
		NodeGammaBeta:               0.2,
		BranchAlpha:                 3,
		BranchBeta:                  100,
		SetSeed:                     19,
		GcContent:                   0.41,
		BiasLeafName:                "",
		NonBiasProbThreshold:        0,
		HighestProbThreshold:        0,
		SeqLen:                      10000,
		SubMatrix:                   true,
	},
	{TestName: "GtrSimGtrRecon",
		SimSubstitutionMatrixFile:   "testdata/substitutionMatrix.txt",
		ReconSubstitutionMatrixFile: "testdata/substitutionMatrix.txt",
		LeavesFile:                  "testdata/leavesFile.txt",
		NewickFile:                  "testdata/currNewick.txt",
		ReconOutFile:                "testdata/reconOutFile.txt",
		NumTrees:                    50,
		NodeGammaAlpha:              3,
		NodeGammaBeta:               0.2,
		BranchAlpha:                 3,
		BranchBeta:                  100,
		SetSeed:                     17,
		GcContent:                   0.41,
		BiasLeafName:                "",
		NonBiasProbThreshold:        0,
		HighestProbThreshold:        0,
		SeqLen:                      10000,
		SubMatrix:                   true,
	},
	{TestName: "GtrSimDefaultRecon",
		SimSubstitutionMatrixFile:   "testdata/substitutionMatrix.txt",
		ReconSubstitutionMatrixFile: "testdata/substitutionMatrix.txt",
		LeavesFile:                  "testdata/leavesFile.txt",
		NewickFile:                  "testdata/currNewick.txt",
		ReconOutFile:                "testdata/reconOutFile.txt",
		NumTrees:                    50,
		NodeGammaAlpha:              3,
		NodeGammaBeta:               0.2,
		BranchAlpha:                 3,
		BranchBeta:                  100,
		SetSeed:                     23,
		GcContent:                   0.41,
		BiasLeafName:                "",
		NonBiasProbThreshold:        0,
		HighestProbThreshold:        0,
		SeqLen:                      10000,
		SubMatrix:                   false,
	},

	{TestName: "JukesSimJukesRecon",
		SimSubstitutionMatrixFile:   "testdata/jukesCantor.txt",
		ReconSubstitutionMatrixFile: "testdata/jukesCantor.txt",
		LeavesFile:                  "testdata/leavesFile.txt",
		NewickFile:                  "testdata/currNewick.txt",
		ReconOutFile:                "testdata/reconOutFile.txt",
		NumTrees:                    50,
		NodeGammaAlpha:              3,
		NodeGammaBeta:               0.2,
		BranchAlpha:                 3,
		BranchBeta:                  100,
		SetSeed:                     29,
		GcContent:                   0.41,
		BiasLeafName:                "",
		NonBiasProbThreshold:        0,
		HighestProbThreshold:        0,
		SeqLen:                      10000,
		SubMatrix:                   true,
	},
	{TestName: "JukesSimGtrRecon",
		SimSubstitutionMatrixFile:   "testdata/jukesCantor.txt",
		ReconSubstitutionMatrixFile: "testdata/substitutionMatrix.txt",
		LeavesFile:                  "testdata/leavesFile.txt",
		NewickFile:                  "testdata/currNewick.txt",
		ReconOutFile:                "testdata/reconOutFile.txt",
		NumTrees:                    50,
		NodeGammaAlpha:              3,
		NodeGammaBeta:               0.2,
		BranchAlpha:                 3,
		BranchBeta:                  100,
		SetSeed:                     31,
		GcContent:                   0.41,
		BiasLeafName:                "",
		NonBiasProbThreshold:        0,
		HighestProbThreshold:        0,
		SeqLen:                      10000,
		SubMatrix:                   true,
	},
	{TestName: "JukesSimDefaultRecon",
		SimSubstitutionMatrixFile:   "testdata/jukesCantor.txt",
		ReconSubstitutionMatrixFile: "testdata/substitutionMatrix.txt",
		LeavesFile:                  "testdata/leavesFile.txt",
		NewickFile:                  "testdata/currNewick.txt",
		ReconOutFile:                "testdata/reconOutFile.txt",
		NumTrees:                    50,
		NodeGammaAlpha:              3,
		NodeGammaBeta:               0.2,
		BranchAlpha:                 3,
		BranchBeta:                  100,
		SetSeed:                     39,
		GcContent:                   0.41,
		BiasLeafName:                "",
		NonBiasProbThreshold:        0,
		HighestProbThreshold:        0,
		SeqLen:                      10000,
		SubMatrix:                   false,
	},
	{TestName: "TransitionSimTransitionRecon",
		SimSubstitutionMatrixFile:   "testdata/extremeSubstitutionMatrix.txt",
		ReconSubstitutionMatrixFile: "testdata/extremeSubstitutionMatrix.txt",
		LeavesFile:                  "testdata/leavesFile.txt",
		NewickFile:                  "testdata/currNewick.txt",
		ReconOutFile:                "testdata/reconOutFile.txt",
		NumTrees:                    50,
		NodeGammaAlpha:              3,
		NodeGammaBeta:               0.2,
		BranchAlpha:                 3,
		BranchBeta:                  100,
		SetSeed:                     39,
		GcContent:                   0.41,
		BiasLeafName:                "",
		NonBiasProbThreshold:        0,
		HighestProbThreshold:        0,
		SeqLen:                      10000,
		SubMatrix:                   true,
	},
}

func TestEmpiricalReconstruction(t *testing.T) {
	var currNumNodes int
	var err error
	var currRandGamma float64
	var currSimTree, currReconTree *expandedTree.ETree
	out := fileio.EasyCreate("testdata/resultsSummary.txt")
	_, err = fmt.Fprintf(out, "Name\tTreeIndex\tNodeName\tInaccuracy\n")
	exception.PanicOnErr(err)
	for _, v := range EmpiricalReconstructionComparison {
		rand.Seed(v.SetSeed)
		for currTreeIndex := 0; currTreeIndex < v.NumTrees; currTreeIndex++ {

			//first, we make a tree to test, and run a molecular evolution simulation to generate sequences
			currRandGamma, _ = numbers.RandGamma(v.NodeGammaAlpha, v.NodeGammaBeta)
			currNumNodes = int(currRandGamma + 2) //ensure we have at least 2 nodes
			if currNumNodes%2 == 0 {              //we want to sample only odd positive integers, this guarantees at lest 3 nodes
				currNumNodes++
			}
			currSimTree = simulate.ETree(currNumNodes, v.BranchAlpha, v.BranchBeta)
			expandedTree.ToNewickFile(v.NewickFile, currSimTree)
			currSimTree.Fasta = &fasta.Fasta{Name: currSimTree.Name, Seq: simulate.RandIntergenicSeq(v.GcContent, v.SeqLen)}
			simulate.NonCoding(currSimTree, v.SimSubstitutionMatrixFile, 0.1)
			WriteLeavesToFasta(currSimTree, v.LeavesFile)

			//second, we run a reconstruction
			currReconTree, err = expandedTree.ReadTree(v.NewickFile, v.LeavesFile)
			exception.PanicOnErr(err)
			unitMatrix := simulate.ParseSubstitutionMatrix(v.ReconSubstitutionMatrixFile)
			expandedTree.PopulateSubstitutionMatrices(currReconTree, unitMatrix, 0.1)
			reconLeaves := expandedTree.GetLeaves(currReconTree)
			for i := range reconLeaves[0].Fasta.Seq {
				LoopNodes(currReconTree, i, v.BiasLeafName, v.NonBiasProbThreshold, v.HighestProbThreshold, v.SubMatrix)
			}

			//third, we compare the reconstruction and sim and write to a file
			currSimBranches := expandedTree.GetBranch(currSimTree)
			reconMap := expandedTree.ToMap(currReconTree)
			for currNode := range currSimBranches {
				_, err = fmt.Fprintf(out, "%v\t%v\t%v\t%v\n", v.TestName, currTreeIndex, currSimBranches[currNode].Name, percentDivergence(currSimBranches[currNode].Fasta.Seq, reconMap[currSimBranches[currNode].Name].Fasta.Seq))
				exception.PanicOnErr(err)
			}
		}
	}
	err = out.Close()
	exception.PanicOnErr(err)
}

func percentDivergence(seqA []dna.Base, seqB []dna.Base) float64 {
	if len(seqA) != len(seqB) {
		log.Fatalf("Error: input sequences are not of the same length.\n")
	}
	var diffCount int = 0
	for currPos := 0; currPos < len(seqA); currPos++ {
		if seqA[currPos] != seqB[currPos] {
			diffCount++
		}
	}
	return float64(diffCount) / float64(len(seqA))
}
*/
