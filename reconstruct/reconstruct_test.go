package reconstruct

import (
	"testing"

	"github.com/vertgenlab/gonomics/dna/pDna"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/expandedTree"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/fasta/pFasta"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/simulate"

	// uncomment for additional tests
	"fmt"
	"log"
	"math/rand"

	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/numbers"
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
	BiasNodeName         string
	NonBiasProbThreshold float64
	HighestProbThreshold float64
	SubMatrix            bool
	PDnaNode             string
	PDnaOutFile          string
	PDnaExpected         string
	Precision            float32
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
		BiasNodeName:         "",
		NonBiasProbThreshold: 0,
		HighestProbThreshold: 0,
		SubMatrix:            false,
		PDnaNode:             "C",
		PDnaOutFile:          "testdata/C.pfa",
		PDnaExpected:         "testdata/CExpected.pfa",
		Precision:            1,
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
	var outPFasta []pFasta.PFasta

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

		if v.PDnaNode != "" {
			outPFasta = []pFasta.PFasta{{Name: v.PDnaNode, Seq: make([]pDna.Float32Base, 0)}}
		}

		for i := 0; i < len(leaves[0].Fasta.Seq); i++ {
			LoopNodes(tree, i, v.BiasLeafName, v.BiasNodeName, v.NonBiasProbThreshold, false, v.HighestProbThreshold, v.SubMatrix, v.PDnaNode, outPFasta)
		}
		WriteTreeToFasta(tree, v.ReconOutFile)

		if v.PDnaNode != "" {
			pFasta.Write(v.PDnaOutFile, outPFasta)
		}

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

		pDnaExpected := pFasta.Read(v.PDnaExpected)
		if !pFasta.AllAreEqual(outPFasta, pDnaExpected, v.Precision) {
			t.Errorf("Error: pFaExtract outFile is not as expected.")
		} else {
			fileio.EasyRemove(v.PDnaOutFile)
		}

		fileio.EasyRemove(v.RandFa)
		fileio.EasyRemove(v.LeavesFile)
		fileio.EasyRemove(v.ReconOutFile)
		fileio.EasyRemove(v.SimTree)
	}
}

// this benchmark takes about 60s CPU to run, so I've commented it out. Go ahead
// and uncomment if you're interested in reproducing it or altering it. -RJM
// /*
// var EmpiricalReconstructionComparison = []struct {
// 	TestName                    string
// 	SimSubstitutionMatrixFile   string
// 	ReconSubstitutionMatrixFile string
// 	LeavesFile                  string
// 	NewickFile                  string
// 	ReconOutFile                string
// 	NumTrees                    int
// 	NodeGammaAlpha              float64
// 	NodeGammaBeta               float64
// 	BranchAlpha                 float64
// 	BranchBeta                  float64
// 	SetSeed                     int64
// 	GcContent                   float64
// 	BiasLeafName                string
// 	BiasN						bool
// 	NonBiasProbThreshold        float64
// 	HighestProbThreshold        float64
// 	SeqLen                      int
// 	SubMatrix                   bool
// 	PDnaNode			 		string
// 	PDnaOutFile			 		string
// 	PDnaExpected		 		string
// }{
// 	{TestName: "GtrSimJukesRecon",
// 		SimSubstitutionMatrixFile:   "testdata/substitutionMatrix.txt",
// 		ReconSubstitutionMatrixFile: "testdata/jukesCantor.txt",
// 		LeavesFile:                  "testdata/leavesFile.txt",
// 		NewickFile:                  "testdata/currNewick.txt",
// 		ReconOutFile:                "testdata/reconOutFile.txt",
// 		NumTrees:                    50,
// 		NodeGammaAlpha:              3,
// 		NodeGammaBeta:               0.2,
// 		BranchAlpha:                 3,
// 		BranchBeta:                  100,
// 		SetSeed:                     19,
// 		GcContent:                   0.41,
// 		BiasLeafName:                "",
// 		BiasN:						 false,
// 		NonBiasProbThreshold:        0,
// 		HighestProbThreshold:        0,
// 		SeqLen:                      10000,
// 		SubMatrix:                   true,
// 	},
// 	{TestName: "GtrSimGtrRecon",
// 		SimSubstitutionMatrixFile:   "testdata/substitutionMatrix.txt",
// 		ReconSubstitutionMatrixFile: "testdata/substitutionMatrix.txt",
// 		LeavesFile:                  "testdata/leavesFile.txt",
// 		NewickFile:                  "testdata/currNewick.txt",
// 		ReconOutFile:                "testdata/reconOutFile.txt",
// 		NumTrees:                    50,
// 		NodeGammaAlpha:              3,
// 		NodeGammaBeta:               0.2,
// 		BranchAlpha:                 3,
// 		BranchBeta:                  100,
// 		SetSeed:                     17,
// 		GcContent:                   0.41,
// 		BiasLeafName:                "",
// 	 	BiasN:						 false,
// 		NonBiasProbThreshold:        0,
// 		HighestProbThreshold:        0,
// 		SeqLen:                      10000,
// 		SubMatrix:                   true,
// 	},
// 	{TestName: "GtrSimDefaultRecon",
// 		SimSubstitutionMatrixFile:   "testdata/substitutionMatrix.txt",
// 		ReconSubstitutionMatrixFile: "testdata/substitutionMatrix.txt",
// 		LeavesFile:                  "testdata/leavesFile.txt",
// 		NewickFile:                  "testdata/currNewick.txt",
// 		ReconOutFile:                "testdata/reconOutFile.txt",
// 		NumTrees:                    50,
// 		NodeGammaAlpha:              3,
// 		NodeGammaBeta:               0.2,
// 		BranchAlpha:                 3,
// 		BranchBeta:                  100,
// 		SetSeed:                     23,
// 		GcContent:                   0.41,
// 		BiasLeafName:                "",
// 	 BiasN:						 false,
// 		NonBiasProbThreshold:        0,
// 		HighestProbThreshold:        0,
// 		SeqLen:                      10000,
// 		SubMatrix:                   false,
// 	},

// 	{TestName: "JukesSimJukesRecon",
// 		SimSubstitutionMatrixFile:   "testdata/jukesCantor.txt",
// 		ReconSubstitutionMatrixFile: "testdata/jukesCantor.txt",
// 		LeavesFile:                  "testdata/leavesFile.txt",
// 		NewickFile:                  "testdata/currNewick.txt",
// 		ReconOutFile:                "testdata/reconOutFile.txt",
// 		NumTrees:                    50,
// 		NodeGammaAlpha:              3,
// 		NodeGammaBeta:               0.2,
// 		BranchAlpha:                 3,
// 		BranchBeta:                  100,
// 		SetSeed:                     29,
// 		GcContent:                   0.41,
// 		BiasLeafName:                "",
// 	 BiasN:						 false,
// 		NonBiasProbThreshold:        0,
// 		HighestProbThreshold:        0,
// 		SeqLen:                      10000,
// 		SubMatrix:                   true,
// 	},
// 	{TestName: "JukesSimGtrRecon",
// 		SimSubstitutionMatrixFile:   "testdata/jukesCantor.txt",
// 		ReconSubstitutionMatrixFile: "testdata/substitutionMatrix.txt",
// 		LeavesFile:                  "testdata/leavesFile.txt",
// 		NewickFile:                  "testdata/currNewick.txt",
// 		ReconOutFile:                "testdata/reconOutFile.txt",
// 		NumTrees:                    50,
// 		NodeGammaAlpha:              3,
// 		NodeGammaBeta:               0.2,
// 		BranchAlpha:                 3,
// 		BranchBeta:                  100,
// 		SetSeed:                     31,
// 		GcContent:                   0.41,
// 		BiasLeafName:                "",
// 	 BiasN:						 false,
// 		NonBiasProbThreshold:        0,
// 		HighestProbThreshold:        0,
// 		SeqLen:                      10000,
// 		SubMatrix:                   true,
// 	},
// 	{TestName: "JukesSimDefaultRecon",
// 		SimSubstitutionMatrixFile:   "testdata/jukesCantor.txt",
// 		ReconSubstitutionMatrixFile: "testdata/substitutionMatrix.txt",
// 		LeavesFile:                  "testdata/leavesFile.txt",
// 		NewickFile:                  "testdata/currNewick.txt",
// 		ReconOutFile:                "testdata/reconOutFile.txt",
// 		NumTrees:                    50,
// 		NodeGammaAlpha:              3,
// 		NodeGammaBeta:               0.2,
// 		BranchAlpha:                 3,
// 		BranchBeta:                  100,
// 		SetSeed:                     39,
// 		GcContent:                   0.41,
// 		BiasLeafName:                "",
// 	 BiasN:						 false,
// 		NonBiasProbThreshold:        0,
// 		HighestProbThreshold:        0,
// 		SeqLen:                      10000,
// 		SubMatrix:                   false,
// 	},
// 	{TestName: "TransitionSimTransitionRecon",
// 		SimSubstitutionMatrixFile:   "testdata/extremeSubstitutionMatrix.txt",
// 		ReconSubstitutionMatrixFile: "testdata/extremeSubstitutionMatrix.txt",
// 		LeavesFile:                  "testdata/leavesFile.txt",
// 		NewickFile:                  "testdata/currNewick.txt",
// 		ReconOutFile:                "testdata/reconOutFile.txt",
// 		NumTrees:                    50,
// 		NodeGammaAlpha:              3,
// 		NodeGammaBeta:               0.2,
// 		BranchAlpha:                 3,
// 		BranchBeta:                  100,
// 		SetSeed:                     39,
// 		GcContent:                   0.41,
// 		BiasLeafName:                "",
// 	 BiasN:						 false,
// 		NonBiasProbThreshold:        0,
// 		HighestProbThreshold:        0,
// 		SeqLen:                      10000,
// 		SubMatrix:                   true,
// 	},
// }

// func TestEmpiricalReconstruction(t *testing.T) {
// 	var currNumNodes int
// 	var err error
// 	var currRandGamma float64
// 	var currSimTree, currReconTree *expandedTree.ETree
// 	var outPFasta []pFasta.PFasta
// 	out := fileio.EasyCreate("testdata/resultsSummary.txt")
// 	_, err = fmt.Fprintf(out, "Name\tTreeIndex\tNodeName\tInaccuracy\n")
// 	exception.PanicOnErr(err)
// 	for _, v := range EmpiricalReconstructionComparison {
// 		rand.New(rand.NewSource(s.SetSeed))
// 		for currTreeIndex := 0; currTreeIndex < v.NumTrees; currTreeIndex++ {

// 			//first, we make a tree to test, and run a molecular evolution simulation to generate sequences
// 			currRandGamma, _ = numbers.RandGamma(v.NodeGammaAlpha, v.NodeGammaBeta)
// 			currNumNodes = int(currRandGamma + 2) //ensure we have at least 2 nodes
// 			if currNumNodes%2 == 0 {              //we want to sample only odd positive integers, this guarantees at lest 3 nodes
// 				currNumNodes++
// 			}
// 			currSimTree = simulate.ETree(currNumNodes, v.BranchAlpha, v.BranchBeta)
// 			expandedTree.ToNewickFile(v.NewickFile, currSimTree)
// 			currSimTree.Fasta = &fasta.Fasta{Name: currSimTree.Name, Seq: simulate.RandIntergenicSeq(v.GcContent, v.SeqLen)}
// 			simulate.NonCoding(currSimTree, v.SimSubstitutionMatrixFile, 0.1)
// 			WriteLeavesToFasta(currSimTree, v.LeavesFile)

// 			//second, we run a reconstruction
// 			currReconTree, err = expandedTree.ReadTree(v.NewickFile, v.LeavesFile)
// 			exception.PanicOnErr(err)
// 			unitMatrix := simulate.ParseSubstitutionMatrix(v.ReconSubstitutionMatrixFile)
// 			expandedTree.PopulateSubstitutionMatrices(currReconTree, unitMatrix, 0.1)
// 			reconLeaves := expandedTree.GetLeaves(currReconTree)
// 			for i := range reconLeaves[0].Fasta.Seq {
// 				LoopNodes(currReconTree, i, v.BiasLeafName, v.NonBiasProbThreshold, v.BiasN, v.HighestProbThreshold, v.SubMatrix, v.PDnaNode, outPFasta)
// 			}

// 			//third, we compare the reconstruction and sim and write to a file
// 			currSimBranches := expandedTree.GetBranch(currSimTree)
// 			reconMap := expandedTree.ToMap(currReconTree)
// 			for currNode := range currSimBranches {
// 				_, err = fmt.Fprintf(out, "%v\t%v\t%v\t%v\n", v.TestName, currTreeIndex, currSimBranches[currNode].Name, percentDivergence(currSimBranches[currNode].Fasta.Seq, reconMap[currSimBranches[currNode].Name].Fasta.Seq))
// 				exception.PanicOnErr(err)
// 			}
// 		}
// 	}
// 	err = out.Close()
// 	exception.PanicOnErr(err)
// }

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

// versions of the above tests, only 1 tree, with pDNA
var EmpiricalReconstructionComparisonPDna = []struct {
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
	BiasNodeName                string
	BiasN                       bool
	NonBiasProbThreshold        float64
	HighestProbThreshold        float64
	SeqLen                      int
	SubMatrix                   bool
	PDnaNode                    string
	PDnaOutFile                 string
	PDnaExpected                string
	Precision                   float32
}{
	{TestName: "GtrSimJukesRecon",
		SimSubstitutionMatrixFile:   "testdata/substitutionMatrix.txt",
		ReconSubstitutionMatrixFile: "testdata/jukesCantor.txt",
		LeavesFile:                  "testdata/leavesFile.txt",
		NewickFile:                  "testdata/currNewick.txt",
		ReconOutFile:                "testdata/reconOutFile.txt",
		NumTrees:                    1,
		NodeGammaAlpha:              3,
		NodeGammaBeta:               0.2,
		BranchAlpha:                 3,
		BranchBeta:                  100,
		SetSeed:                     19,
		GcContent:                   0.41,
		BiasLeafName:                "",
		BiasN:                       false,
		NonBiasProbThreshold:        0,
		HighestProbThreshold:        0,
		SeqLen:                      10000,
		SubMatrix:                   true,
		PDnaNode:                    "Child_4",
		PDnaOutFile:                 "testdata/GtrSimJukesReconChild4.pfa",
		PDnaExpected:                "testdata/GtrSimJukesExpectedChild4.pfa",
		Precision:                   1e-3,
	},
	{TestName: "GtrSimGtrRecon",
		SimSubstitutionMatrixFile:   "testdata/substitutionMatrix.txt",
		ReconSubstitutionMatrixFile: "testdata/substitutionMatrix.txt",
		LeavesFile:                  "testdata/leavesFile.txt",
		NewickFile:                  "testdata/currNewick.txt",
		ReconOutFile:                "testdata/reconOutFile.txt",
		NumTrees:                    1,
		NodeGammaAlpha:              3,
		NodeGammaBeta:               0.2,
		BranchAlpha:                 3,
		BranchBeta:                  100,
		SetSeed:                     17,
		GcContent:                   0.41,
		BiasLeafName:                "",
		BiasN:                       false,
		NonBiasProbThreshold:        0,
		HighestProbThreshold:        0,
		SeqLen:                      10000,
		SubMatrix:                   true,
		PDnaNode:                    "Child_3",
		PDnaOutFile:                 "testdata/GtrSimGtrReconChild3.pfa",
		PDnaExpected:                "testdata/GtrSimGtrExpectedChild3.pfa",
		Precision:                   1e-3,
	},
	{TestName: "GtrSimDefaultRecon",
		SimSubstitutionMatrixFile:   "testdata/substitutionMatrix.txt",
		ReconSubstitutionMatrixFile: "testdata/substitutionMatrix.txt",
		LeavesFile:                  "testdata/leavesFile.txt",
		NewickFile:                  "testdata/currNewick.txt",
		ReconOutFile:                "testdata/reconOutFile.txt",
		NumTrees:                    1,
		NodeGammaAlpha:              3,
		NodeGammaBeta:               0.2,
		BranchAlpha:                 3,
		BranchBeta:                  100,
		SetSeed:                     23,
		GcContent:                   0.41,
		BiasLeafName:                "",
		BiasN:                       false,
		NonBiasProbThreshold:        0,
		HighestProbThreshold:        0,
		SeqLen:                      10000,
		SubMatrix:                   false,
		PDnaNode:                    "Child_10",
		PDnaOutFile:                 "testdata/GtrSimDefaultReconChild10.pfa",
		PDnaExpected:                "testdata/GtrSimDefaultExpectedChild10.pfa",
		Precision:                   1e-3,
	},
	{TestName: "JukesSimJukesRecon",
		SimSubstitutionMatrixFile:   "testdata/jukesCantor.txt",
		ReconSubstitutionMatrixFile: "testdata/jukesCantor.txt",
		LeavesFile:                  "testdata/leavesFile.txt",
		NewickFile:                  "testdata/currNewick.txt",
		ReconOutFile:                "testdata/reconOutFile.txt",
		NumTrees:                    1,
		NodeGammaAlpha:              3,
		NodeGammaBeta:               0.2,
		BranchAlpha:                 3,
		BranchBeta:                  100,
		SetSeed:                     29,
		GcContent:                   0.41,
		BiasLeafName:                "",
		BiasN:                       false,
		NonBiasProbThreshold:        0,
		HighestProbThreshold:        0,
		SeqLen:                      10000,
		SubMatrix:                   true,
		PDnaNode:                    "Child_16",
		PDnaOutFile:                 "testdata/JukesSimJukesReconChild16.pfa",
		PDnaExpected:                "testdata/JukesSimJukesExpectedChild16.pfa",
		Precision:                   1e-3,
	},
	{TestName: "JukesSimGtrRecon",
		SimSubstitutionMatrixFile:   "testdata/jukesCantor.txt",
		ReconSubstitutionMatrixFile: "testdata/substitutionMatrix.txt",
		LeavesFile:                  "testdata/leavesFile.txt",
		NewickFile:                  "testdata/currNewick.txt",
		ReconOutFile:                "testdata/reconOutFile.txt",
		NumTrees:                    1,
		NodeGammaAlpha:              3,
		NodeGammaBeta:               0.2,
		BranchAlpha:                 3,
		BranchBeta:                  100,
		SetSeed:                     31,
		GcContent:                   0.41,
		BiasLeafName:                "",
		BiasN:                       false,
		NonBiasProbThreshold:        0,
		HighestProbThreshold:        0,
		SeqLen:                      10000,
		SubMatrix:                   true,
		PDnaNode:                    "Child_25",
		PDnaOutFile:                 "testdata/JukesSimGtrReconChild25.pfa",
		PDnaExpected:                "testdata/JukesSimGtrExpectedChild25.pfa",
		Precision:                   1e-3,
	},
	{TestName: "JukesSimDefaultRecon",
		SimSubstitutionMatrixFile:   "testdata/jukesCantor.txt",
		ReconSubstitutionMatrixFile: "testdata/substitutionMatrix.txt",
		LeavesFile:                  "testdata/leavesFile.txt",
		NewickFile:                  "testdata/currNewick.txt",
		ReconOutFile:                "testdata/reconOutFile.txt",
		NumTrees:                    1,
		NodeGammaAlpha:              3,
		NodeGammaBeta:               0.2,
		BranchAlpha:                 3,
		BranchBeta:                  100,
		SetSeed:                     39,
		GcContent:                   0.41,
		BiasLeafName:                "",
		BiasN:                       false,
		NonBiasProbThreshold:        0,
		HighestProbThreshold:        0,
		SeqLen:                      10000,
		SubMatrix:                   false,
		PDnaNode:                    "Child_4",
		PDnaOutFile:                 "testdata/JukesSimDefaultReconChild4.pfa",
		PDnaExpected:                "testdata/JukesSimDefaultExpectedChild4.pfa",
		Precision:                   1e-3,
	},
	{TestName: "TransitionSimTransitionRecon",
		SimSubstitutionMatrixFile:   "testdata/extremeSubstitutionMatrix.txt",
		ReconSubstitutionMatrixFile: "testdata/extremeSubstitutionMatrix.txt",
		LeavesFile:                  "testdata/leavesFile.txt",
		NewickFile:                  "testdata/currNewick.txt",
		ReconOutFile:                "testdata/reconOutFile.txt",
		NumTrees:                    1,
		NodeGammaAlpha:              3,
		NodeGammaBeta:               0.2,
		BranchAlpha:                 3,
		BranchBeta:                  100,
		SetSeed:                     39,
		GcContent:                   0.41,
		BiasLeafName:                "",
		BiasN:                       false,
		NonBiasProbThreshold:        0,
		HighestProbThreshold:        0,
		SeqLen:                      10000,
		SubMatrix:                   true,
		PDnaNode:                    "Child_5",
		PDnaOutFile:                 "testdata/TransitionSimTransitionReconChild5.pfa",
		PDnaExpected:                "testdata/TransitionSimTransitionExpectedChild5.pfa",
		Precision:                   1e-3,
	},
}

func TestEmpiricalReconstruction(t *testing.T) {
	var currNumNodes int
	var err error
	var currRandGamma float64
	var currSimTree, currReconTree *expandedTree.ETree
	var outPFasta []pFasta.PFasta
	out := fileio.EasyCreate("testdata/resultsSummary.txt")
	_, err = fmt.Fprintf(out, "Name\tTreeIndex\tNodeName\tInaccuracy\n")
	exception.PanicOnErr(err)
	for _, v := range EmpiricalReconstructionComparisonPDna {
		rand.New(rand.NewSource(v.SetSeed))

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

		// initialise output pFasta
		if v.PDnaNode != "" {
			outPFasta = []pFasta.PFasta{{Name: v.PDnaNode, Seq: make([]pDna.Float32Base, 0)}}
		}

		//second, we run a reconstruction
		currReconTree, err = expandedTree.ReadTree(v.NewickFile, v.LeavesFile)
		exception.PanicOnErr(err)
		unitMatrix := simulate.ParseSubstitutionMatrix(v.ReconSubstitutionMatrixFile)
		expandedTree.PopulateSubstitutionMatrices(currReconTree, unitMatrix, 0.1)
		reconLeaves := expandedTree.GetLeaves(currReconTree)
		for i := range reconLeaves[0].Fasta.Seq {
			LoopNodes(currReconTree, i, v.BiasLeafName, v.BiasNodeName, v.NonBiasProbThreshold, v.BiasN, v.HighestProbThreshold, v.SubMatrix, v.PDnaNode, outPFasta)
		}

		//third, we compare the reconstruction and sim and write to a file
		currSimBranches := expandedTree.GetBranch(currSimTree)
		reconMap := expandedTree.ToMap(currReconTree)
		for currNode := range currSimBranches {
			_, err = fmt.Fprintf(out, "%v\t%v\t%v\n", v.TestName, currSimBranches[currNode].Name, percentDivergence(currSimBranches[currNode].Fasta.Seq, reconMap[currSimBranches[currNode].Name].Fasta.Seq))
			exception.PanicOnErr(err)
		}

		// fourth, compare the reconstructed pfasta and write to file
		if v.PDnaNode != "" {
			pFasta.Write(v.PDnaOutFile, outPFasta)
			expectedPfa := pFasta.Read(v.PDnaExpected)
			if !pFasta.AllAreEqual(expectedPfa, outPFasta, v.Precision) {
				fmt.Fprintf(out, "%v\t%v\n", v.TestName, v.PDnaNode)
			} else {
				fileio.EasyRemove(v.PDnaOutFile)
			}
		}
	}
	err = out.Close()
	exception.PanicOnErr(err)
}
