package main

import (
	"flag"
	"fmt"
	"log"
	"math/rand"
	"os"

	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/expandedTree"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/numbers"
	"github.com/vertgenlab/gonomics/simulate"
)

// NonCodingSettings defines usage settings for the simulateEvol nonCoding subcommand.
type NonCodingSettings struct {
	TreeFile               string
	FastaFile              string
	OutFile                string
	UnitBranchLength       float64
	SubstitutionMatrixFile string
	NumNodes               int
	GammaAlpha             float64
	GammaBeta              float64
	GcContent              float64
	LenSeq                 int
	SetSeed                int64
	NewickOut              string
}

// NonCodingUsage defines the usage statement for the simulateEvol nonCoding subcommand.
func NonCodingUsage(nonCodingFlags *flag.FlagSet) {
	fmt.Print(
		"simulateEvol nonCoding - a tool to simulate molecular evolution in noncoding regions.\n" +
			"This program simulates Jukes-Cantor evolution by default, but accepts custom substitution matrices.\n" +
			"This program does not support indels, but rather simulates substitutions.\n" +
			"The program can simulate molecular evolution along a specified input Newick tree, or can generate its own.\n" +
			"Similarly, a root DNA sequence can be provided in fasta format, or a root sequence can be simulated.\n" +
			"Usage:\n" +
			"\tsimulateEvol nonCoding out.fasta\n" +
			"options:\n")
	nonCodingFlags.PrintDefaults()
}

// parseNonCodingArgs is the main function of the simulateEvol nonCoding subcommand. It parses options and launches the NonCoding function.
func parseNonCodingArgs() {
	var expectedNumArgs int = 1
	var err error
	nonCodingFlags := flag.NewFlagSet("nonCoding", flag.ExitOnError)
	nonCodingFlags.Usage = func() { NonCodingUsage(nonCodingFlags) }
	var setSeed *int64 = nonCodingFlags.Int64("setSeed", -1, "Use a specific seed for the RNG.")
	var unitBranchLength *float64 = nonCodingFlags.Float64("unitBranchLength", -100, "Set the branch length over which a custom substitution matrix was derived. If ")
	var substitutionMatrix *string = nonCodingFlags.String("substitutionMatrixFile", "", "Specify a custom substitution matrix.")
	var numNodes *int = nonCodingFlags.Int("numNodes", 13, "If generating a Newick tree, set the total number of nodes in the simulate tree.")
	var gammaAlpha *float64 = nonCodingFlags.Float64("gammaAlpha", 1, "If generating a Newick tree, set the alpha parameter for Gamma-distributed branch lengths.")
	var gammaBeta *float64 = nonCodingFlags.Float64("gammaBeta", 50, "If generating a Newick tree, set the beta parameter for Gamma-distributed branch lengths.")
	var gcContent *float64 = nonCodingFlags.Float64("gcContent", 0.41, "If generating a root DNA sequence, set the GC content for the simulated sequence.")
	var lenSeq *int = nonCodingFlags.Int("lenSeq", 100, "If generating a root DNA sequence, set the length of the simulated sequence.")
	var treeFile *string = nonCodingFlags.String("treeFile", "", "Specify a file for simulating molecular evolution along a pre-specified Newick tree.")
	var fastaFile *string = nonCodingFlags.String("fastaFile", "", "Specify a sequence for the root node of the simulation. This file is expected to contain only one sequence. Output name is hardcoded as 'root', so the fasta name will be ignored.")
	var newickOut *string = nonCodingFlags.String("newickOut", "", "Write the tree to an output Newick-format file.")
	err = nonCodingFlags.Parse(os.Args[2:])
	exception.PanicOnErr(err)
	if len(nonCodingFlags.Args()) != expectedNumArgs {
		nonCodingFlags.Usage()
		log.Fatalf("Error: expecting %d arguments, but got %d\n",
			expectedNumArgs, len(nonCodingFlags.Args()))
	}
	outFile := nonCodingFlags.Arg(0)
	s := NonCodingSettings{
		TreeFile:               *treeFile,
		FastaFile:              *fastaFile,
		OutFile:                outFile,
		SetSeed:                *setSeed,
		NumNodes:               *numNodes,
		GammaAlpha:             *gammaAlpha,
		GammaBeta:              *gammaBeta,
		GcContent:              *gcContent,
		LenSeq:                 *lenSeq,
		SubstitutionMatrixFile: *substitutionMatrix,
		UnitBranchLength:       *unitBranchLength,
		NewickOut:              *newickOut,
	}
	NonCoding(s)
}

// NonCoding simulates molecular evolution along a Newick tree and writes the resulting sequences to a file.
// A Newick tree file can be provided. Alternatively, one can be generated with a user-specified number of
// nodes and Gamma-distribute random branch lengths.
func NonCoding(s NonCodingSettings) {
	rand.New(rand.NewSource(s.SetSeed))
	var answer []fasta.Fasta
	var root *expandedTree.ETree
	var err error

	if s.GammaAlpha <= 0 || s.GammaBeta <= 0 {
		log.Fatalf("Error: expected Gamma distribution parameters to be positive numbers. Found: alpha=%v, beta=%v\n", s.GammaAlpha, s.GammaBeta)
	}
	if s.GcContent < 0 || s.GcContent > 1 {
		log.Fatalf("Error: GcContent must be a value between 0 and 1. Found: %v.\n", s.GcContent)
	}
	if s.LenSeq < 0 {
		log.Fatalf("Error: expected lenSeq to be a positive number. Found: %v.\n", s.LenSeq)
	}
	if s.TreeFile != "" {
		root, err = expandedTree.ReadNewick(s.TreeFile)
	} else {
		root = simulate.ETree(s.NumNodes, s.GammaAlpha, s.GammaBeta)
	}
	if s.UnitBranchLength < 0 {
		s.UnitBranchLength, _ = numbers.RandGamma(s.GammaAlpha, s.GammaBeta)
	}
	if s.FastaFile != "" {
		records := fasta.Read(s.FastaFile)
		if len(records) != 1 {
			log.Fatalf("Error: expected 1 sequence in the input fasta file. Received: %v.\n", len(records))
		}
		fasta.ToUpper(records[0])
		root.Fasta = &records[0]
		root.Name = "root"
	} else {
		root.Fasta = &fasta.Fasta{Name: "root", Seq: simulate.RandIntergenicSeq(s.GcContent, s.LenSeq)}
	}
	exception.PanicOnErr(err)
	root = simulate.NonCoding(root, s.SubstitutionMatrixFile, s.UnitBranchLength)
	nodes := expandedTree.GetTree(root)
	for currNode := range nodes {
		answer = append(answer, *nodes[currNode].Fasta)
	}
	fasta.Write(s.OutFile, answer)
	if s.NewickOut != "" {
		out := fileio.EasyCreate(s.NewickOut)
		_, err = fmt.Fprintf(out, "%s\n", expandedTree.ToNewickString(root))
		exception.PanicOnErr(err)
		err = out.Close()
		exception.PanicOnErr(err)
	}
}
