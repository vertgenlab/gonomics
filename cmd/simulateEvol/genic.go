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
	"github.com/vertgenlab/gonomics/simulate"
)

type GenicSettings struct {
	TreeFile     string
	InFile       string
	OutFile      string
	SetSeed      int64
	SimOutFile   string
	GenePredFile string
}

func GenicUsage(genicFlags *flag.FlagSet) {
	fmt.Print(
		"simulateEvol genic - a tool to simulate molecular evolution in coding regions, using the BLOSOM matrix.\n" +
			"Usage:\n" +
			"simulateEvol genic tree.newick in.fasta out.fasta\n" +
			"options:\n")
	genicFlags.PrintDefaults()
}

func parseGenicArgs() {
	var expectedNumArgs int = 3
	var err error
	genicFlags := flag.NewFlagSet("genic", flag.ExitOnError)
	genicFlags.Usage = func() { GenicUsage(genicFlags) }
	var setSeed *int64 = genicFlags.Int64("setSeed", -1, "Use a specific seed for the RNG.")
	var completeSimOutputFile *string = genicFlags.String("simOutFile", "", "Also return a fasta file containing the sequence for every node in the tree, including internal nodes.")
	var genePredFile *string = genicFlags.String("genePredFile", "", "Specify gene features. Sequence evolution in gene features will follow the BLOSOM matrix.")
	err = genicFlags.Parse(os.Args[2:])
	exception.PanicOnErr(err)
	if len(genicFlags.Args()) != expectedNumArgs {
		genicFlags.Usage()
		log.Fatalf("Error: expecting %d arguments, but got %d\n",
			expectedNumArgs, len(genicFlags.Args()))
	}

	treeFile := genicFlags.Arg(0)
	inFile := genicFlags.Arg(1)
	outFile := genicFlags.Arg(2)

	s := GenicSettings{
		TreeFile:     treeFile,
		InFile:       inFile,
		OutFile:      outFile,
		SetSeed:      *setSeed,
		SimOutFile:   *completeSimOutputFile,
		GenePredFile: *genePredFile,
	}
	Genic(s)
}

func Genic(s GenicSettings) {
	var fastas, leafFastas []fasta.Fasta
	rng := rand.New(rand.NewSource(s.SetSeed))

	tree, err := expandedTree.ReadTree(s.TreeFile, s.InFile)
	exception.PanicOnErr(err)
	simulate.Simulate(s.InFile, tree, s.GenePredFile, true, rng)
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
	fasta.Write(s.OutFile, leafFastas)
}
