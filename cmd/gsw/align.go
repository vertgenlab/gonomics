package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/align"
	"github.com/vertgenlab/gonomics/chromInfo"
	"github.com/vertgenlab/gonomics/sam"
	"github.com/vertgenlab/gonomics/simpleGraph"
	"os"
	"strings"
)

type GswSettings struct {
	Cmd      *flag.FlagSet
	Index    int
	StepSize int
	Cpus     int
	Matrix   string
	Liftover string
	Help     string
	Out      string
}

func alignUsage() {
	fmt.Printf(
		"  align\t\tGraph-Smith-Waterman: align single or paired end fastqs\n")
}

func extendedAlignUsage() {
	fmt.Print(
		"Usage:\n" +
			"  gsw align [options] ref.gg R1.fastq.gz R2.fastq.gz\n\n" +
			"Required:\n" +
			"  ref[.gg/.fa]\t\tReference genome: both fasta and graph formats accepted\n\n" +
			"Options:\n" +
			"  -i, --index\t\tKmer length for index hash look up  (default: 32)\n" +
			"  -s, --step\t\tOffset sliding window step size for hash set up (default 32)\n" +
			"  -t, --threads\t\tNumber of CPUs for Goroutine concurrency (default: 4)\n" +
			"  -p, --project\t\tConvert alignment coordinate to linear reference in sam format\n" +
			"  -m, --matrix\t\tScores used to align matches and mismatches (default: humanChimp)\n\n")
}

func initArgsGsw() *GswSettings {
	gsw := &GswSettings{Cmd: flag.NewFlagSet("align", flag.ExitOnError)}
	gsw.Cmd.Usage = extendedAlignUsage

	gsw.Cmd.IntVar(&gsw.Index, "index", 32, "``hash look up length")
	gsw.Cmd.IntVar(&gsw.Index, "i", 32, "``hash look up length")

	gsw.Cmd.IntVar(&gsw.StepSize, "window", 32, "offset position of sliding window of hash")
	gsw.Cmd.IntVar(&gsw.StepSize, "w", 32, "offset position of sliding window of hash")

	gsw.Cmd.IntVar(&gsw.Cpus, "threads", 4, "Number of CPUs for goroutines")
	gsw.Cmd.IntVar(&gsw.Cpus, "t", 4, "Number of CPUs for goroutines")

	gsw.Cmd.StringVar(&gsw.Matrix, "matrix", "humanChimp", "Scoring matrix for alignment")
	gsw.Cmd.StringVar(&gsw.Matrix, "m", "humanChimp", "Scoring matrix for alignment")

	gsw.Cmd.StringVar(&gsw.Liftover, "liftover", "", "liftover to linear reference sam file")
	gsw.Cmd.StringVar(&gsw.Liftover, "l", "", "liftover to linear reference sam file")

	gsw.Cmd.StringVar(&gsw.Help, "help", "", "display usage")
	gsw.Cmd.StringVar(&gsw.Help, "h", "", "display usage")

	gsw.Cmd.StringVar(&gsw.Out, "out", "/dev/stdout", "Output filename, [.gg/.vcf/.sam]")
	gsw.Cmd.StringVar(&gsw.Out, "o", "/dev/stdout", "Output filename, [.gg/.vcf/.sam]")

	return gsw
}

func RunAlignExe() {
	gsw := initArgsGsw()
	gsw.Cmd.Parse(os.Args[2:])
	if len(gsw.Cmd.Args()) == 0 {
		gsw.Cmd.Usage()
	} else {
		graphSmithWaterman(gsw.Index, gsw.StepSize, gsw.Cpus, gsw.Matrix, gsw.Out, gsw.Liftover, gsw.Cmd.Args())
	}
}

//TODO: add check for too many args or an implementation for multiple fastq pairs
func graphSmithWaterman(seedNum int, stepSize int, cpus int, score string, out string, liftover string, args []string) {
	//should be at most 3 args to add to the input, reference, readOne and/or readTwo
	var genomeGraph *simpleGraph.SimpleGraph = simpleGraph.Read(args[0])
	var readOne string
	var readTwo string

	if len(args) > 1 {
		readOne = args[1]
	}
	if len(args) == 3 {
		readTwo = args[2]
	}
	switch true {
	case len(args) == 2 && !strings.HasSuffix(liftover, ".sizes"):
		GswToGiraf(genomeGraph, readOne, out, cpus, seedNum, stepSize, selectScoreMatrix(score))
	case len(args) == 3 && !strings.HasSuffix(liftover, ".sizes"):
		GswToGirafPair(genomeGraph, readOne, readTwo, out, cpus, seedNum, stepSize, selectScoreMatrix(score))
	case strings.HasSuffix(liftover, ".sizes"):
		chrSize := chromInfo.ReadToSlice(liftover)
		header := sam.ChromInfoSamHeader(chrSize)
		if len(args) == 2 {
			GswToSam(genomeGraph, readOne, out, cpus, seedNum, stepSize, selectScoreMatrix(score), header)
		}
		if len(args) == 3 {
			GswToSamPair(genomeGraph, readOne, readTwo, out, cpus, seedNum, stepSize, selectScoreMatrix(score), header)
		}
	default:
		extendedAlignUsage()
		errorMessage()
	}
}

func selectScoreMatrix(score string) [][]int64 {
	var scoreMatrix [][]int64
	switch {
	case strings.Contains(score, "humanChimp"):
		scoreMatrix = align.HumanChimpTwoScoreMatrix
	case strings.Contains(score, "hoxD55"):
		scoreMatrix = align.HoxD55ScoreMatrix
	case strings.Contains(score, "mouseRat"):
		scoreMatrix = align.MouseRatScoreMatrix
	case strings.Contains(score, "general"):
		scoreMatrix = align.DefaultScoreMatrix
	default:
		errorMessage()
	}
	return scoreMatrix
}
