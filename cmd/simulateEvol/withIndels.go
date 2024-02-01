package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/simulate"
	"log"
	"math/rand"
	"os"
)

type WithIndelsSettings struct {
	FastaFile      string
	OutFile        string
	BranchLength   float64
	PropIndels     float64
	Lambda         float64
	GcContent      float64
	TransitionBias float64
	VcfOutFile     string
	QName          string
	SetSeed        int64
}

func WithIndelsUsage(withIndelsFlags *flag.FlagSet) {
	fmt.Print(
		"simulateEvol withIndels - simulates molecular evolution with\n" +
			"a two-parameter mutation model and geometrically-distributed INDELs.\n" +
			"Usage:\n" +
			"simulateEvol withIndels seq.fasta outFile.fasta\n" +
			"options:\n")
	withIndelsFlags.PrintDefaults()
}

func parseWithIndelsFlags() {
	var expectedNumArgs int = 2
	var err error
	withIndelsFlags := flag.NewFlagSet("withIndels", flag.ExitOnError)
	withIndelsFlags.Usage = func() { WithIndelsUsage(withIndelsFlags) }
	var lambda *float64 = withIndelsFlags.Float64("lambda", 1, "Define the rate parameter for the exponential distribution used to generate simulated INDEL lengths.")
	var propIndel *float64 = withIndelsFlags.Float64("propIndel", 0, "Proportion of simulated variants that should be insertions or deletions.")
	var branchLength *float64 = withIndelsFlags.Float64("branchLength", 0, "Specify the divergence rate (must be between 0 and 1).")
	var gcContent *float64 = withIndelsFlags.Float64("gcContent", 0.42, "Set the GC content for simulated insertion sequences.")
	var setSeed *int64 = withIndelsFlags.Int64("setSeed", -1, "Use a specific seed for the RNG.")

	err = withIndelsFlags.Parse(os.Args[2:])
	exception.PanicOnErr(err)
	withIndelsFlags.Usage = func() { WithIndelsUsage(withIndelsFlags) }

	if len(withIndelsFlags.Args()) != expectedNumArgs {
		withIndelsFlags.Usage()
		log.Fatalf("Error: expecting %d arguments, but got %d\n",
			expectedNumArgs, len(withIndelsFlags.Args()))
	}

	inFile := withIndelsFlags.Arg(0)
	outFile := withIndelsFlags.Arg(1)

	s := WithIndelsSettings{
		FastaFile:    inFile,
		OutFile:      outFile,
		Lambda:       *lambda,
		PropIndels:   *propIndel,
		BranchLength: *branchLength,
		GcContent:    *gcContent,
		SetSeed:      *setSeed,
	}
	WithIndels(s)
}

func WithIndels(s WithIndelsSettings) {
	if s.PropIndels < 0 || s.PropIndels > 1 {
		log.Fatalf("The propIndels option must be a value between 0 and 1.")
	}
	if s.GcContent < 0 || s.GcContent > 1 {
		log.Fatalf("GcContent must be a value between 0 and 1.")
	}
	if s.TransitionBias < 0 {
		log.Fatalf("TransitionBias must be a nonnegative number.")
	}
	if s.BranchLength < 0 || s.BranchLength > 1 {
		log.Fatalf("The branchLength argument must be a value between 0 and 1.")
	}
	rand.Seed(s.SetSeed)
	outFasta := simulate.WithIndels(s.FastaFile, s.BranchLength, s.PropIndels, s.Lambda, s.GcContent, s.TransitionBias, s.VcfOutFile, s.QName)
	fasta.Write(s.OutFile, outFasta)
}
