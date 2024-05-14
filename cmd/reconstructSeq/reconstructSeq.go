// Command Group: "Sequence Evolution & Reconstruction"

// Reconstruct ancient sequences using extant genomes and a newick tree with branch lengths
package main

import (
	"flag"
	"fmt"
	"log"
	"strings"

	"github.com/vertgenlab/gonomics/dna/pDna"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/expandedTree"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/fasta/pFasta"
	"github.com/vertgenlab/gonomics/reconstruct"
	"github.com/vertgenlab/gonomics/simulate"
)

type Settings struct {
	NewickInput            string
	FastaInput             string
	OutFile                string
	BiasLeafName           string
	BiasNodeName           string
	NonBiasProbThreshold   float64
	BiasN                  bool
	HighestProbThreshold   float64
	KeepAllSeq             bool
	SubMatrix              bool
	UnitBranchLength       float64
	SubstitutionMatrixFile string
	PDnaNode               string
	PDnaNodeMulti          []string
	PDnaOutFile            string
}

func ReconstructSeq(s Settings) {
	var treeFastas []fasta.Fasta

	if s.NonBiasProbThreshold < 0 || s.NonBiasProbThreshold > 1 {
		log.Fatalf("Error: nonBiasProbThreshold must be a value between 0 and 1. Found: %v.\n", s.NonBiasProbThreshold)
	}

	if s.NonBiasProbThreshold > 0 && s.BiasLeafName == "" {
		log.Fatalf("Error: nonBiasProbThreshold was set, but no BiasLeafName was provided.\n")
	}

	if s.HighestProbThreshold < 0 || s.HighestProbThreshold > 1 {
		log.Fatalf("Error: highestProbThreshold must be a value between 0 and 1. Found: %v.\n", s.HighestProbThreshold)
	}

	tree, err := expandedTree.ReadTree(s.NewickInput, s.FastaInput)
	exception.FatalOnErr(err)
	if s.SubMatrix {
		unitMatrix := simulate.ParseSubstitutionMatrix(s.SubstitutionMatrixFile)
		expandedTree.PopulateSubstitutionMatrices(tree, unitMatrix, s.UnitBranchLength)
	}

	leaves := expandedTree.GetLeaves(tree)
	branches := expandedTree.GetBranch(tree)

	//initialise pfasta for pdna node
	pDnaRecords := []pFasta.PFasta{{Name: s.PDnaNode, Seq: make([]pDna.Float32Base, 0)}}
	// use len(s.PDnaNodeMulti) > 0 to check for PDnaNodeMulti mode and initialize multi pfasta if necessary
	var pDnaRecordsMulti []pFasta.PFasta
	if len(s.PDnaNodeMulti) > 0 {
		// have reference sequence at the top of output pfasta
		refPfa := pFasta.FaToPfa(*leaves[0].Fasta) // convert records[0] aka leaves[0].Fasta from Fasta to Pfasta
		pDnaRecordsMulti = append(pDnaRecordsMulti, refPfa)
		// start pDnaNodeMulti entries in pDnaRecordsMulti, add names now, append sequences in reconstruct.LoopNodes
		for _, v := range s.PDnaNodeMulti {
			pDnaRecordsMulti = append(pDnaRecordsMulti, pFasta.PFasta{Name: v})
		}
	}

	for i := range leaves[0].Fasta.Seq {
		// TODO: need to fix LoopNodes calls in other packages
		reconstruct.LoopNodes(tree, i, s.BiasLeafName, s.BiasNodeName, s.NonBiasProbThreshold, s.BiasN, s.HighestProbThreshold, s.SubMatrix, s.PDnaNode, s.PDnaNodeMulti, pDnaRecords, pDnaRecordsMulti)
	}

	for j := range leaves {
		treeFastas = append(treeFastas, *leaves[j].Fasta)
	}

	for k := range branches {
		treeFastas = append(treeFastas, *branches[k].Fasta)
	}

	if s.KeepAllSeq { // in the keepAllSeq option
		records := fasta.Read(s.FastaInput) // fasta.Read already makes sure that sequence names are unique
		treeFastasMap := fasta.ToMap(treeFastas)
		var found bool
		for i := range records {
			_, found = treeFastasMap[records[i].Name]
			if !found {
				if i == 0 { // if reference fasta is not in tree, then make it the first fasta, before the treeFastas
					treeFastas = append([]fasta.Fasta{records[i]}, treeFastas...) // this can be achieved by making records[i] a slice, and appending "treeFastas..."
				} else {
					treeFastas = append(treeFastas, records[i]) // if non-reference fasta is not in tree, then append it after the treeFastas
				}
			}
		}
	}

	fasta.Write(s.OutFile, treeFastas)

	if s.PDnaNode != "" {
		pFasta.Write(s.PDnaOutFile, pDnaRecords)
	}
	if len(s.PDnaNodeMulti) > 0 {
		pFasta.Write(s.PDnaOutFile, pDnaRecordsMulti)
	}
}

func usage() {
	fmt.Print(
		"reconstructSeq performs ancestral sequence reconstruction based on an input multiFa alignment and Newick tree." +
			"This program returns a fasta file containing sequences for all nodes of the tree, including the input sequences (leaves)," +
			"and the inferred ancestral nodes.\n" +
			"Usage:\n" +
			"reconstructSeq <treeFile.txt> <inFasta.fasta> <outFasta.fasta>\n" +
			"options:\n")
	flag.PrintDefaults()
}

func main() {
	var biasLeafName *string = flag.String("biasLeafName", "", "Specify an extant (leaf) sequence towards which we will bias reconstruction for the immediate ancestor (parent node), unless another node is specified.")
	var biasNodeName *string = flag.String("biasNodeName", "", "Specify the node for biased reconstruction.")
	var nonBiasProbThreshold *float64 = flag.Float64("nonBiasBaseThreshold", 0, "Given that a biasLeafName specifies a reference species, when reconstructing the sequence of a non-reference species, unless the sum of probabilities for all non-reference bases is above this value, the reference base is returned.")
	var biasN *bool = flag.Bool("biasN", false, "Given -biasLeafName and -nonBiasBasThreshold, unless the sum of probabilities for all non-reference bases is above this value, dna.N is returned.")
	var highestProbThreshold *float64 = flag.Float64("highestProbThreshold", 0, "The highest probability base must be above this value to be accepted for reconstruction. Otherwise, dna.N will be returned.")
	var keepAllSeq *bool = flag.Bool("keepAllSeq", false, "By default, reconstructSeq discards sequences in the fasta input that are not specified in the newick input, because they are not used in the reconstruction. If keepAllSeq is set to TRUE, reconstructSeq will keep all sequences in the fasta input, even if they are not used in the reconstruction.")
	var substitutionMatrixFile *string = flag.String("substitutionMatrixFile", "", "Set a file to define a substitution matrix.")
	var unitBranchLength *float64 = flag.Float64("unitBranchLength", -1, "If using a substitution matrix, specify the branch length over which the substitution matrix was derived.")
	var subMatrix *bool = flag.Bool("subMatrix", false, "Use a substitution matrix instead of the default model. If no substitution matrix file is provided, the Jukes-Cantor model will be used.")
	var pDnaNode *string = flag.String("pDnaNode", "", "Specify a node to get the pDNA sequence of. Defaults to empty. Requires pDnaOutFile")
	var pDnaNodeMulti *string = flag.String("pDnaNodeMulti", "", "Specify >1 nodes in a comma-delimited list to get the pDNA sequence of. e.g. -pDnaNodeMulti=hca,hoa. The input multiFa alignment's reference sequence (aka sequence index 0, first sequence) will be converted from fasta to pFasta and become the reference sequence in the output pFasta. Defaults to empty. Requires pDnaOutFile")
	var pDnaOutFile *string = flag.String("pDnaOutFile", "", "Name of pDnaNode pfasta. Requires pDnaNode or pDnaNodeMulti.")

	var expectedNumArgs = 3

	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime | log.Lshortfile)
	flag.Parse()

	if len(flag.Args()) != expectedNumArgs {
		flag.Usage()
		log.Fatalf("Error: expecting %d arguments, but got %d\n",
			expectedNumArgs, len(flag.Args()))
	}

	var pDnaNodeMultiSplit []string
	if *pDnaNodeMulti != "" {
		pDnaNodeMultiSplit = strings.Split(*pDnaNodeMulti, ",")
	}

	newickInput := flag.Arg(0)
	fastaInput := flag.Arg(1)
	outFile := flag.Arg(2)

	s := Settings{
		NewickInput:            newickInput,
		FastaInput:             fastaInput,
		OutFile:                outFile,
		BiasLeafName:           *biasLeafName,
		BiasNodeName:           *biasNodeName,
		NonBiasProbThreshold:   *nonBiasProbThreshold,
		BiasN:                  *biasN,
		HighestProbThreshold:   *highestProbThreshold,
		KeepAllSeq:             *keepAllSeq,
		SubstitutionMatrixFile: *substitutionMatrixFile,
		UnitBranchLength:       *unitBranchLength,
		SubMatrix:              *subMatrix,
		PDnaNode:               *pDnaNode,
		PDnaNodeMulti:          pDnaNodeMultiSplit,
		PDnaOutFile:            *pDnaOutFile,
	}

	ReconstructSeq(s)
}
