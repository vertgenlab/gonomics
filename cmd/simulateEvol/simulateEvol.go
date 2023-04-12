// Command Group: "Sequence Evolution & Reconstruction"

package main

import (
	"flag"
	"fmt"
	"log"
	"math/rand"

	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/expandedTree"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/simulate"
)

type Settings struct {
	FastaFile      string
	TreeFile       string
	LeafOutFile    string
	GenePredFile   string
	SimOutFile     string
	Lambda         float64
	PropIndel      float64
	BranchLength   float64
	SetSeed        int64
	GcContent      float64
	VcfOutFile     string
	TransitionBias float64
	QName          string
}

func SimulateEvol(s Settings) {
	rand.Seed(s.SetSeed)
	var fastas, leafFastas []fasta.Fasta

	if s.BranchLength < 0 || s.BranchLength > 1 {
		log.Fatalf("The branchLength argument must be a value between 0 and 1.")
	}
	if s.BranchLength != 0 && s.TreeFile != "" {
		log.Fatalf("The user must choose either branchLength mode or newick mode, both cannot be selected at once.")
	}
	if s.BranchLength == 0 && s.TreeFile == "" {
		log.Fatalf("The user must choose either branchLength mode or newick mode, neither option has been selected.")
	}
	if s.PropIndel < 0 || s.PropIndel > 1 {
		log.Fatalf("The propIndel option must be a value between 0 and 1.")
	}
	if s.GcContent < 0 || s.GcContent > 1 {
		log.Fatalf("GcContent must be a value between 0 and 1.")
	}
	if s.TransitionBias < 0 {
		log.Fatalf("TransitionBias must be a nonnegative number.")
	}

	if s.TreeFile != "" {
		tree, err := expandedTree.ReadTree(s.TreeFile, s.FastaFile)
		exception.PanicOnErr(err)
		simulate.Simulate(s.FastaFile, tree, s.GenePredFile, true)
		nodes := expandedTree.GetTree(tree)

		for i := 0; i < len(nodes); i++ {
			fastas = append(fastas, *nodes[i].Fasta)
			if nodes[i].Left == nil && nodes[i].Right == nil {
				leafFastas = append(leafFastas, *nodes[i].Fasta)
			}
		}
		if s.SimOutFile != "" {
			fasta.Write(s.SimOutFile, fastas)
		}
	} else {
		leafFastas = simulate.SimulateWithIndels(s.FastaFile, s.BranchLength, s.PropIndel, s.Lambda, s.GcContent, s.TransitionBias, s.VcfOutFile, s.QName)
	}
	fasta.Write(s.LeafOutFile, leafFastas)

}

func usage() {
	fmt.Print(
		"simulateEvol simulates sequence evolution on an input fasta.\n" +
			"The user can provide one of two options.\n" +
			"First, a user may use the branchLength option to generate a single simulated sequence from the input fasta.\n" +
			"Second, a user may use the newickTree option to simulate sequences across the tree, treating the input sequence as the root node.\n" +
			"Note that the newickTree mode does not support insertions, but may introduce deletions." +
			"By default, newickTree mode returns only the leaf sequences in fasta format.\n" +
			"Usage:\n" +
			"simulateEvol <rootFasta.fasta> <outputFile.fasta> \n" +
			"options:\n")
	flag.PrintDefaults()
}

func main() {
	var expectedNumArgs = 2
	var newickTree *string = flag.String("newickTree", "", "Specify a newick tree filename to simulate evolution for each node.")
	var completeSimOutputFile *string = flag.String("simOutFile", "", "Also return a fasta file containing the sequence for every node in the tree, including internal nodes.")
	var genePredFile *string = flag.String("genePredFile", "", "Specify gene features. Sequence evolution in gene features will follow the BLOSOM matrix.")
	var lambda *float64 = flag.Float64("lambda", 1, "Define the rate parameter for the exponential distribution used to generate simulated INDEL lengths.")
	var propIndel *float64 = flag.Float64("propIndel", 0, "Proportion of simulated variants that should be insertions or deletions.")
	var branchLength *float64 = flag.Float64("branchLength", 0, "For a single descendent sequence, specify the divergence rate (must be between 0 and 1).")
	var gcContent *float64 = flag.Float64("gcContent", 0.42, "Set the GC content for simulated insertion sequences.")
	var setSeed *int64 = flag.Int64("setSeed", -1, "Use a specific seed for the RNG.")
	var vcfOutFile *string = flag.String("vcfOutFile", "", "For branchLength mode, specify a vcf filename to record simulated mutations.")
	var transitionBias *float64 = flag.Float64("transitionBias", 1, "Set a bias for transitions over transversions during sequence evolution. Defaults to the Jukes-Cantor model, where the transition bias is 1 (even with transversion frequency).")
	var qName *string = flag.String("qName", "evol", "Set the suffix for the output sequence fasta name. Default suffix evol for an example chr1 will appear as chr1_evol.")

	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime | log.Lshortfile)
	flag.Parse()

	if len(flag.Args()) != expectedNumArgs {
		flag.Usage()
		log.Fatalf("Error: expecting %d arguments, but got %d\n",
			expectedNumArgs, len(flag.Args()))
	}

	rootFasta := flag.Arg(0)
	leafOutFile := flag.Arg(1)

	s := Settings{
		FastaFile:      rootFasta,
		TreeFile:       *newickTree,
		LeafOutFile:    leafOutFile,
		SimOutFile:     *completeSimOutputFile,
		GenePredFile:   *genePredFile,
		BranchLength:   *branchLength,
		Lambda:         *lambda,
		PropIndel:      *propIndel,
		GcContent:      *gcContent,
		SetSeed:        *setSeed,
		VcfOutFile:     *vcfOutFile,
		TransitionBias: *transitionBias,
		QName:          *qName,
	}

	SimulateEvol(s)
}
