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

// this part is from pfaFindFast.go
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

	// Conversion to valid sequences. TODO: remove after debugging?
	firstQuery = pFasta.MakeValid(firstQuery)
	secondQuery = pFasta.MakeValid(secondQuery)

	// QC firstQuery and secondQuery sequences. TODO: remove after debugging?
	//pFasta.QC(firstQuery)
	//pFasta.QC(secondQuery)
	fmt.Printf("QC complete\n")
	//fmt.Printf("firstQuery: %v\n", firstQuery)
	//fmt.Printf("secondQuery: %v\n", secondQuery)
	speedyWindowDifference(reference, firstQuery, secondQuery, s)
}

// part from pfaFindFast.go ended

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
	// TODO: a better way may be to separate ReconstructSeq and ReconstructSeqPfa into sub-functions?
	// min requires go 1.21
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

	// keepAllSeq for pFasta too
	//fmt.Printf("I think this is reference sequence: leaves[0].Fasta.Name: %v\n", leaves[0].Fasta.Name) // confirmed humanT2T is reference sequence
	refPfa := pFasta.FaToPfa(*leaves[0].Fasta) // convert records[0] aka leaves[0].Fasta from Fasta to Pfasta? Is this necessary?
	answerPfa = append([]pFasta.PFasta{refPfa}, answerPfa...)

	fasta.Write(s.OutFile, treeFastas)

	if s.OutPfaFile != "" {
		// bypass writing/reading pFasta file. TODO: better to fix writing/reading and have separate reconstructSeq and pFafindfast cmd
		s2.InPfasta = answerPfa
		pfaFindFast(s2)

		for _, v := range answerPfa { // instead of printing whole pfa, which creates large file that can't be read properly, just QC each base
			for i := range v.Seq {
				// TODO: check invalid base before writing pFastsa. Remove after debugging
				if math.IsNaN(float64(v.Seq[i].A)) || math.IsNaN(float64(v.Seq[i].C)) || math.IsNaN(float64(v.Seq[i].G)) || math.IsNaN(float64(v.Seq[i].T)) || v.Seq[i].A < 0 || v.Seq[i].C < 0 || v.Seq[i].G < 0 || v.Seq[i].T < 0 {
					log.Fatalf("Before writing pFasta and found invalid base: %v at position %v, and fasta name is: %v\n", v.Seq[i], i, v.Name)
				}
				// no if print all on short example
				if i == 10 || i == 5000 || i == 5102 || i == 5120 {
					fmt.Printf("Before writing pFasta, Check base: %v at position %v, and fasta name is: %v\n", v.Seq[i], i, v.Name)
				}
			}
			//fmt.Printf("%v: %v\n", i, v)
		}

		pFasta.Write(s.OutPfaFile, answerPfa)
		//TODO: put the below as another function in pfa package to print human-readable (non-binary) pfa to terminal or remove/comment-out after debugging
		outpfa := pFasta.Read(s.OutPfaFile)
		for _, v := range outpfa { // instead of printing whole pfa, which creates large file that can't be read properly, just QC each base
			for i := range v.Seq {
				// TODO: check invalid base when reading pFastsa. Remove after debugging
				if math.IsNaN(float64(v.Seq[i].A)) || math.IsNaN(float64(v.Seq[i].C)) || math.IsNaN(float64(v.Seq[i].G)) || math.IsNaN(float64(v.Seq[i].T)) || v.Seq[i].A < 0 || v.Seq[i].C < 0 || v.Seq[i].G < 0 || v.Seq[i].T < 0 {
					log.Fatalf("Read pFasta and found invalid base: %v at position %v, and fasta name is: %v\n", v.Seq[i], i, v.Name)
				}
				// no if print all on short example
				if i == 10 || i == 5000 || i == 5102 || i == 5120 {
					fmt.Printf("Read pFasta, Check base: %v at position %v, and fasta name is: %v\n", v.Seq[i], i, v.Name)
				}
			}
			//fmt.Printf("%v: %v\n", i, v)
		}
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
	var outPfaFile *string = flag.String("outPfaFile", "", "Write a pFasta file for the inferred ancestral nodes specified by name in the -pfaNames option.")
	var pfaNames *string = flag.String("pfaNames", "", "A comma-delimited list of the names of inferred ancestral nodes to write a pFasta file for, e.g. -pfaNames=hca,hoa.")
	// this part is from pfaFindFast.go
	var pfafindfastOutfile *string = flag.String("pfafindfastOutfile", "pfafindfastOutfile.bed", "Specify output of pfafindfast.")
	var firstQueryName *string = flag.String("firstQueryName", "", "Specify the name of the first query sequence")
	var secondQueryName *string = flag.String("secondQueryName", "", "Specify the name of the second query sequence")
	var windowSize *int = flag.Int("windowSize", 1000, "Specify the window size")
	var refChromName *string = flag.String("chrom", "", "Specify a chrom name of the reference sequence")
	var removeN *bool = flag.Bool("removeN", false, "Excludes bed regions with Ns in the reference from the output.")
	var longOutput *bool = flag.Bool("longOutput", false, "Print percent diverged and raw -Log10PValue in output. Requires the 'divergenceRate' argument.")
	var divergenceRate *float64 = flag.Float64("divergenceRate", math.MaxFloat64, "Set the null divergence rate for p value calculations with 'longOutput'.")
	var outputAlnPos *bool = flag.Bool("outputAlnPos", false, "Print the alignment position of the window's start in output as the last column.")
	// part from pfaFindFast.go ended

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

	// part from pfafindfast.go
	if *longOutput && *divergenceRate == math.MaxFloat64 {
		log.Fatalf("Error: must set a 'divergenceRate' if using the 'longOutput' option.\n")
	}

	if *divergenceRate != math.MaxFloat64 {
		if *divergenceRate < 0 || *divergenceRate > 1 {
			log.Fatalf("Error: divergence rate must be a value between 0 and 1.\n")
		}
	}
	// part from pfafindfast.go ended

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

	// part from pfafindfast.go
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
	// part from pfafindfast.go ended

	ReconstructSeq(s, s2)
}
