package main

import(
	"flag"
	"os"
	"fmt"
	"strings"
	"github.com/vertgenlab/gonomics/sam"
	"github.com/vertgenlab/gonomics/simpleGraph"
	"github.com/vertgenlab/gonomics/chromInfo"
	"github.com/vertgenlab/gonomics/align"
)

type AlignExe struct {
	Cmd *flag.FlagSet
	Index int
	StepSize int
	Cpus int
	Matrix string
	Liftover string
	Help string
	Out string
}

func Gsw() *AlignExe {

	gsw := &AlignExe{Cmd: flag.NewFlagSet("align", flag.ContinueOnError)}

	gsw.Cmd.IntVar(&gsw.Index, "index", 32, "``hash look up length")
	gsw.Cmd.IntVar(&gsw.Index, "i", 32, "``hash look up length")

	gsw.Cmd.IntVar(&gsw.StepSize, "window", 32, "offset position of sliding window of hash")
	gsw.Cmd.IntVar(&gsw.StepSize, "w", 32, "offset position of sliding window of hash")
	
	gsw.Cmd.IntVar(&gsw.Cpus, "threads", 4, "Number of CPUs for goroutines")
	gsw.Cmd.IntVar(&gsw.Cpus, "t", 4, "Number of CPUs for goroutines")
	
	gsw.Cmd.StringVar(&gsw.Matrix,"matrix", "humanChimp", "Scoring matrix for alignment")
	gsw.Cmd.StringVar(&gsw.Matrix,"m", "humanChimp", "Scoring matrix for alignment")

	gsw.Cmd.StringVar(&gsw.Liftover, "liftover", "", "liftover to linear reference sam file")
	gsw.Cmd.StringVar(&gsw.Liftover, "l", "", "liftover to linear reference sam file")

	gsw.Cmd.StringVar(&gsw.Help, "help", "", "display usage")
	gsw.Cmd.StringVar(&gsw.Help, "h", "", "display usage")

	gsw.Cmd.StringVar(&gsw.Out, "out", "/dev/stdout", "Output filename, [.gg/.vcf/.sam]")
	gsw.Cmd.StringVar(&gsw.Out, "o", "/dev/stdout", "Output filename, [.gg/.vcf/.sam]")
	
	gsw.Cmd.Usage = alignExtend
	return gsw
}

func RunAlignExe() error {
	
	gsw := Gsw()

	Init(gsw.Cmd, os.Args[2:])
	tail := gsw.Cmd.Args()
	if len(tail) == 0 || gsw.Help != "" {
		flag.PrintDefaults()
		return fmt.Errorf("Error: Apologies, your command prompt was not recognized...\n\n-xoxo GG\n")
	} else {
		graphSmithWaterman(gsw.Index, gsw.StepSize, gsw.Cpus, gsw.Matrix, gsw.Out, gsw.Liftover, tail)
		return nil
	}
}

func graphSmithWaterman(seedNum int, stepSize int, cpus int, score string, out string, liftover string, tail []string) {
	gRef := simpleGraph.Read(tail[0])
	switch len(tail) {
	case 2:
		if strings.HasSuffix(liftover, ".sizes") {
			chrSize := chromInfo.ReadToSlice(liftover)
			header := sam.ChromInfoSamHeader(chrSize)
			WrapSingleGirafLiftover(gRef, tail[1], out, cpus, seedNum, stepSize, selectScoreMatrix(score), header)
		} else {
			GswToGiraf(gRef, tail[1], out, cpus, seedNum, stepSize, selectScoreMatrix(score))
		}
	case 3:
		if strings.HasSuffix(liftover, ".sizes") {
			chrSize := chromInfo.ReadToSlice(liftover)
			header := sam.ChromInfoSamHeader(chrSize)
			WrapGirafLiftoverToSam(gRef, tail[1], tail[2], out, cpus, seedNum, stepSize, selectScoreMatrix(score), header)
		} else {
			GswToGirafPair(gRef, tail[1], tail[2], out, cpus, seedNum, stepSize, selectScoreMatrix(score))
		}
	default:
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
	}
	return scoreMatrix
}
