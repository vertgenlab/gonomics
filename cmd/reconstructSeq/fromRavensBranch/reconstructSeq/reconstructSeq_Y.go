// Command Group: "Sequence Evolution & Reconstruction"

// Reconstruct ancient sequences using extant genomes and a newick tree with branch lengths
package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/dna/pDna"
	"github.com/vertgenlab/gonomics/fasta/pFasta"
	"github.com/vertgenlab/gonomics/simulate"
	"log"
	"math"
	"strings"

	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/expandedTree"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/reconstruct"
)

// this type block is from pfaFindFast.go
type Settings2 struct {
	InPfasta        []pFasta.PFasta
	OutFile         string
	FirstQueryName  string
	SecondQueryName string
	WindowSize      int
	RefChromName    string
	RemoveN         bool
	LongOutput      bool
	DivergenceRate  float64
	OutputAlnPos    bool
}

// this func block is from pfaFindFast.go
func pfaFindFast(s Settings2) {
	var reference, firstQuery, secondQuery []pDna.Float32Base
	var found bool

	records := s.InPfasta               // fasta.Read already makes sure that sequence names are unique
	recordsMap := pFasta.ToMap(records) // generate map[string][]dna.Base. No order but can do key-value search

	if len(records) < 2 {
		log.Fatalf("Error: There must be at least 2 pFasta records in the input file.\n")
	}

	//if reference and query names were not specified in the command, then use 1st and 2nd fasta records' names, respectively
	//allows backward compatibility with previously-written code
	if s.FirstQueryName == "" {
		firstQuery = records[0].Seq
	} else {
		_, found = recordsMap[s.FirstQueryName]
		if found {
			firstQuery = recordsMap[s.FirstQueryName]
		} else {
			log.Fatalf("Error: first query name is not found in the input file.\n")
		}
	}
	if s.SecondQueryName == "" {
		secondQuery = records[1].Seq
	} else {
		_, found = recordsMap[s.SecondQueryName]
		if found {
			secondQuery = recordsMap[s.SecondQueryName]
		} else {
			log.Fatalf("Error: second query name is not found in the input file.\n")
		}
	}

	reference = records[0].Seq // reference will always be the 1st pFasta record in a multi-pFa input file

	if !(len(reference) == len(firstQuery) && len(reference) == len(secondQuery)) {
		log.Fatalf("Error: Reference, first query, and second query sequences are not all of equal length.\n")
	}

	speedyWindowDifference(reference, firstQuery, secondQuery, s)
}

type Settings struct {
	NewickInput            string
	FastaInput             string
	OutFile                string
	BiasLeafName           string
	NonBiasProbThreshold   float64
	HighestProbThreshold   float64
	KeepAllSeq             bool
	SubMatrix              bool
	UnitBranchLength       float64
	SubstitutionMatrixFile string
	OutPfaFile             string
	PfaNames               []string
}

func ReconstructSeq(s Settings, s2 Settings2) {
	var treeFastas []fasta.Fasta
	var answerPfa = make([]pFasta.PFasta, max(len(s.PfaNames), 1))

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

	for i := range leaves[0].Fasta.Seq {
		if s.OutPfaFile != "" {
			answerPfa = reconstruct.LoopNodesPfa(tree, i, s.BiasLeafName, s.NonBiasProbThreshold, s.HighestProbThreshold, answerPfa, s.PfaNames)
		} else {
			reconstruct.LoopNodes(tree, i, s.BiasLeafName, s.NonBiasProbThreshold, s.HighestProbThreshold, s.SubMatrix)
		}
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

	if s.OutPfaFile != "" {
		// must have reference as Pfasta first sequence
		refPfa := pFasta.FaToPfa(*leaves[0].Fasta) // convert records[0] aka leaves[0].Fasta from Fasta to Pfasta
		answerPfa = append([]pFasta.PFasta{refPfa}, answerPfa...)
	}

	fasta.Write(s.OutFile, treeFastas)

	if s.OutPfaFile != "" {
		pFasta.Write(s.OutPfaFile, answerPfa)

		// these 2 lines bypass writing and reading pFasta file, directly pipe reconstructSeq into pfaFindFast
		s2.InPfasta = answerPfa
		pfaFindFast(s2)
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
	var biasLeafName *string = flag.String("biasLeafName", "", "Specify an extant (leaf) sequence towards which we will bias reconstruction for the immediate ancestor (parent node).")
	var nonBiasProbThreshold *float64 = flag.Float64("nonBiasBaseThreshold", 0, "Given that a biasLeafName specifies a reference species, when reconstructing the sequence of a non-reference species, unless the sum of probabilities for all non-reference bases is above this value, the reference base is returned.")
	var highestProbThreshold *float64 = flag.Float64("highestProbThreshold", 0, "The highest probability base must be above this value to be accepted for reconstruction. Otherwise, dna.N will be returned.")
	var keepAllSeq *bool = flag.Bool("keepAllSeq", false, "By default, reconstructSeq discards sequences in the fasta input that are not specified in the newick input, because they are not used in the reconstruction. If keepAllSeq is set to TRUE, reconstructSeq will keep all sequences in the fasta input, even if they are not used in the reconstruction.")
	var substitutionMatrixFile *string = flag.String("substitutionMatrixFile", "", "Set a file to define a substitution matrix.")
	var unitBranchLength *float64 = flag.Float64("unitBranchLength", -1, "If using a substitution matrix, specify the branch length over which the substitution matrix was derived.")
	var subMatrix *bool = flag.Bool("subMatrix", false, "Use a substitution matrix instead of the default model. If no substitution matrix file is provided, the Jukes-Cantor model will be used.")
	var outPfaFile *string = flag.String("outPfaFile", "", "Write one pFasta file for the inferred ancestral nodes specified by name in the -pfaNames option. Input multiFa alignment's reference sequence (aka sequence index 0, first sequence) will be converted from fasta to pFasta and become the reference sequence in the output pFasta.")
	var pfaNames *string = flag.String("pfaNames", "", "A comma-delimited list of the names of inferred ancestral nodes to write a pFasta file for, e.g. -pfaNames=hca,hoa. Input multiFa alignment's reference sequence (aka sequence index 0, first sequence) will be converted from fasta to pFasta and become the reference sequence in the output pFasta.")
	// this part (below lines of var declaration) is from pfaFindFast.go
	var pfafindfastOutfile *string = flag.String("pfafindfastOutfile", "pfafindfastOutfile.bed", "Specify output of pfafindfast.")
	var firstQueryName *string = flag.String("firstQueryName", "", "Specify the name of the first query sequence")
	var secondQueryName *string = flag.String("secondQueryName", "", "Specify the name of the second query sequence")
	var windowSize *int = flag.Int("windowSize", 1000, "Specify the window size")
	var refChromName *string = flag.String("chrom", "", "Specify a chrom name of the reference sequence")
	var removeN *bool = flag.Bool("removeN", false, "Excludes bed regions with Ns in the reference from the output.")
	var longOutput *bool = flag.Bool("longOutput", false, "Print percent diverged and raw -Log10PValue in output. Requires the 'divergenceRate' argument.")
	var divergenceRate *float64 = flag.Float64("divergenceRate", math.MaxFloat64, "Set the null divergence rate for p value calculations with 'longOutput'.")
	var outputAlnPos *bool = flag.Bool("outputAlnPos", false, "Print the alignment position of the window's start in output as the last column.")

	var expectedNumArgs = 3

	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime | log.Lshortfile)
	flag.Parse()

	if len(flag.Args()) != expectedNumArgs {
		flag.Usage()
		log.Fatalf("Error: expecting %d arguments, but got %d\n",
			expectedNumArgs, len(flag.Args()))
	}

	var pfaNamesSplit []string
	if *outPfaFile != "" {
		if *pfaNames == "" {
			flag.Usage()
			log.Fatalf("Error: -outPfaFile was specified, but -pfaNames was empty\n")
		} else { // -outPfaFile and -pfaNames were both specified. Parse -pfaNames
			pfaNamesSplit = strings.Split(*pfaNames, ",")
		}
	} else if *outPfaFile == "" && *pfaNames != "" {
		flag.Usage()
		log.Fatalf("Error: -outPfaFile was not specified, but -pfaNames was not empty\n")
	} // else is -outPfaFile and -pfaNames were both empty, aka mode without pFa

	// this if check is from pfaFindFast.go
	if *longOutput && *divergenceRate == math.MaxFloat64 {
		log.Fatalf("Error: must set a 'divergenceRate' if using the 'longOutput' option.\n")
	}

	// this if check is from pfaFindFast.go
	if *divergenceRate != math.MaxFloat64 {
		if *divergenceRate < 0 || *divergenceRate > 1 {
			log.Fatalf("Error: divergence rate must be a value between 0 and 1.\n")
		}
	}

	newickInput := flag.Arg(0)
	fastaInput := flag.Arg(1)
	outFile := flag.Arg(2)

	s := Settings{
		NewickInput:            newickInput,
		FastaInput:             fastaInput,
		OutFile:                outFile,
		BiasLeafName:           *biasLeafName,
		NonBiasProbThreshold:   *nonBiasProbThreshold,
		HighestProbThreshold:   *highestProbThreshold,
		KeepAllSeq:             *keepAllSeq,
		SubstitutionMatrixFile: *substitutionMatrixFile,
		UnitBranchLength:       *unitBranchLength,
		SubMatrix:              *subMatrix,
		OutPfaFile:             *outPfaFile,
		PfaNames:               pfaNamesSplit,
	}

	// this var declaration block is from pfaFindFast.go
	s2 := Settings2{
		InPfasta:        []pFasta.PFasta{},
		OutFile:         *pfafindfastOutfile,
		FirstQueryName:  *firstQueryName,
		SecondQueryName: *secondQueryName,
		WindowSize:      *windowSize,
		RefChromName:    *refChromName,
		RemoveN:         *removeN,
		LongOutput:      *longOutput,
		DivergenceRate:  *divergenceRate,
		OutputAlnPos:    *outputAlnPos,
	}

	ReconstructSeq(s, s2)
}
