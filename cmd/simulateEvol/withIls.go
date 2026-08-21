package main

import (
	"flag"
	"fmt"
	"log"
	"os"

	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/expandedTree"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/numbers/matrix"
	"github.com/vertgenlab/gonomics/simulate"
)

// IlsSettings defines usage settings for the simulateEvol ils subcommand.
type IlsSettings struct {
	RootsFile              string
	TransitionMatrixFile   string
	ChromName              string
	OutPathPrefix          string
	UnitBranchLength       float64
	AncSeqFile             string
	LenSeq                 int64
	SetSeed                int64
	LeafFastasOnly         bool
	SubstitutionMatrixFile string
}

// IlsUsage defines the usage statement for the simulateEvol ils subcommand.
func IlsUsage(ilsFlags *flag.FlagSet) {
	fmt.Print(
		"simulateEvol ils - simulate evolution with incomplete lineage sorting.\n" +
			"This program simulates molecular evolution along a specified set of input Newick trees." +
			"This program simulates with Jukes-Cantor evolution by default, but accepts custom substitution matrices.\n" +
			"This program does not support indels, but rather simulates substitutions.\n" +
			"The program can take in a specified ancestral sequence or randomly generate an initial ancestral sequence\n" +
			"Usage:\n" +
			"\tsimulateEvol ils roots.txt transition_matrix.tsv anc.fasta outPathPrefix unitBranchLength \n" +
			"options:\n",
	)
	ilsFlags.PrintDefaults()
}

// transitionmatrix, chromName, outpathprefix, and unit branch length are mandatory.
// leaffastas only if not provided is assumed false.
// either anc seq or lenseq + setseed must be provided

// parseIlsArgs is the main function of the simulateEvol nonCoding subcommand. It parses options and launches the NonCoding function.
func parseIlsArgs() {
	var expectedNumArgs int = 1
	var err error
	ilsFlags := flag.NewFlagSet("ils", flag.ExitOnError)
	ilsFlags.Usage = func() { NonCodingUsage(ilsFlags) }

	var rootsFile *string = ilsFlags.String("rootsFile", "", "Specify a file for simulating molecular evolution along a set of pre-specified Newick trees.")
	var transitionMatrixFile *string = ilsFlags.String("transitionMatrixFile", "", "Specify a file that describes the probability of transitions between topology states.")
	var ancSeqFile *string = ilsFlags.String("ancSeqFile", "", "Specify the initial ancestral sequence. If empty, must provide setSEed and lenSeq.")
	var setSeed *int64 = ilsFlags.Int64("setSeed", -1, "Use a specific seed for the RNG.")
	var lenSeq *int64 = ilsFlags.Int64("lenSeq", -1, "If generating a root DNA sequence, set the length of the simulated sequence. Ignored if ancSeqFile provided.")
	var chromName *string = ilsFlags.String("chromName", "", "Specify the name of the output sequence.")
	var outPathPrefix *string = ilsFlags.String("outPathPrefix", "", "Specify the output directory and prefix of output files.")
	var leafFastasOnly *bool = ilsFlags.Bool("outPathPrefix", false, "Specify if only leaf fastas are provided in output. Defaults to false.")
	var substitutionMatrixFile *string = ilsFlags.String("substitutionMatrixFile", "", "Specify a custom substitution matrix.")
	var unitBranchLength *float64 = ilsFlags.Float64("unitBranchLength", -1, "Set the branch length over which a custom substitution matrix was derived.")

	err = ilsFlags.Parse(os.Args[2:])
	exception.PanicOnErr(err)
	if len(ilsFlags.Args()) != expectedNumArgs {
		ilsFlags.Usage()
		log.Fatalf("Error: expecting %d arguments, but got %d\n",
			expectedNumArgs, len(ilsFlags.Args()))
	}

	s := IlsSettings{
		RootsFile:              *rootsFile,
		TransitionMatrixFile:   *transitionMatrixFile,
		ChromName:              *chromName,
		OutPathPrefix:          *outPathPrefix,
		UnitBranchLength:       *unitBranchLength,
		AncSeqFile:             *ancSeqFile,
		LenSeq:                 *lenSeq,
		SetSeed:                *setSeed,
		LeafFastasOnly:         *leafFastasOnly,
		SubstitutionMatrixFile: *substitutionMatrixFile,
	}
	Ils(s)
}

// Ils simulates molecular evolution along a Newick tree and writes the resulting sequences to a file.
// A Newick tree file can be provided. Alternatively, one can be generated with a user-specified number of
// nodes and Gamma-distribute random branch lengths.
func Ils(s IlsSettings) {
	var err error

	m, err := matrix.ReadDense(s.TransitionMatrixFile, '\t')
	if err != nil {
		log.Fatalf("error reading transition matrix: %v", err)
	}

	if len(s.RootsFile) == 0 {
		log.Fatalf("Must provide roots file.")
	} else if len(s.TransitionMatrixFile) == 0 {
		log.Fatalf("Must provide transition file.")
	} else if len(s.ChromName) == 0 {
		log.Fatalf("Must provide desired output sequence name.")
	} else if s.UnitBranchLength == -1 {
		log.Fatalf("Must provide unit branch length as used in roots files.")
	}

	var ancSeq []fasta.Fasta
	if s.AncSeqFile != "" {
		ancSeq = fasta.Read(s.AncSeqFile)
	} else {
		if s.LenSeq == -1 {
			log.Fatalf("Must provide either ancestral sequence or desired length of randomly generated sequence.")
		}
		ancSeq = []fasta.Fasta{}
	}

	// you need to read in the s.RootsFile and then make that rootsFiles = []string{}
	rootsFiles := fileio.Read(s.RootsFile)
	roots := make([]*expandedTree.ETree, len(rootsFiles))
	for i, rootFile := range rootsFiles {
		var err error
		roots[i], err = expandedTree.ReadNewick(rootFile)
		if err != nil {
			log.Fatalf("Error: could not read Newick file at %q: %v\n", rootFile, err)
		}
		exception.PanicOnErr(err)
	}

	if s.SetSeed == -1 {
		log.Fatalf("Must provide random seed.")
	}

	anc, evolved, topoRecord, ilsEvolved := simulate.SimulateIls(roots, m, ancSeq, int(s.LenSeq), s.SetSeed, s.ChromName, s.LeafFastasOnly, s.SubstitutionMatrixFile, s.UnitBranchLength)

	fasta.Write(fmt.Sprintf("%s_anc.fasta", s.OutPathPrefix), anc)
	for idx, rec := range evolved {
		fasta.Write(fmt.Sprintf("%s_forward_evolved_topo_v%d.fasta", s.OutPathPrefix, idx), rec)
	}
	bed.Write(fmt.Sprintf("%s.bed", s.OutPathPrefix), topoRecord)
	fasta.Write(fmt.Sprintf("%s_ils.fasta", s.OutPathPrefix), ilsEvolved)
}
