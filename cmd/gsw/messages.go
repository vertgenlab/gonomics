package main

import (
	"fmt"
	//"github.com/vertgenlab/gonomics/align"
	"github.com/vertgenlab/gonomics/sam"
	"github.com/vertgenlab/gonomics/simpleGraph"
	"strings"
	"flag"
	"log"
	"os"
)
//Genome Graph Gonomics
func usage() {
	fmt.Print(
		"GSW - Graph Smith Waterman:\n\nGonomics Genome Graph Software\n" +
			"Author: Eric Au \teric.au@duke.edu\n\tCraig Lowe\tcraig.lowe@duke.edu\n\n" +
			"Source code: https://github.com/vertgenlab/gonomics\n" +
			"Documents: https://github.com/edotau/sticklebackCipher\n\n" +
			"Version: 0.1.0\n\n" +
			"Usage:\n" +
			"  gsw [options]\n\n" +
			"Options:\n")
	alignUsage()
	ggToolsUsage()
	viewUsage()
	helpMessgae()
	flagsPrint()
}

var extendHelpMsg *flag.FlagSet = flag.NewFlagSet("help", flag.ExitOnError)

func helpMessgae() {
	fmt.Printf(
		"  help\t\tDetailed help message for any command\n")
}
func alignUsage() {
	fmt.Printf(
		"  align\t\tGraph-Smith-Waterman: align single or paired end fastqs\n")
}
func ggToolsUsage() {
	fmt.Printf(
		"  ggtools\tGenomic utilities to create, manipulate and operate on genome graphs\n")
}
func viewUsage() {
	fmt.Printf(
		"  view\t\tVisualize graph generated alignment\n")
}
func flagsPrint() {
	fmt.Print(
		"\nFlags:\n" +
			"  -h, --help\t\tEnter gsw --help [align/ggtools/view] for detailed information\n" +
			"  -o, --out\t\tFilename[.gg/.vcf/.gz/.sam]  (default: /dev/stdout)\n\n")
	//"  -t, --threads\t\tNumber of CPUs for goroutines  (default: 4)\n\n")
}

func errorMessage() {
	usage()
	//fmt.Print("\n\n")
	log.Fatalf("Error: Apologies, your command prompt was not recognized...\n\n-xoxo GG\n")
}

func moreHelp(cmdFlag string) {
	if strings.Contains("align", cmdFlag) {
		alignExtend()
	} else if strings.Contains("ggtools", cmdFlag) {
		ggtoolsExtend()
	} else if strings.Contains("view", cmdFlag){
		viewExtend()
	} else {
		errorMessage()
	}
}

func alignExtend() {
	fmt.Print(
		"Usage:\n" +
			"  gsw align [options] ref.gg R1.fastq.gz R2.fastq.gz\n\n" +
			"Required:\n" +
			"  ref[.gg/.fa]\t\tReference genome: acceptes both fasta and graph formats\n\n" +
			"Options:\n" +
			"  -i, --index\t\tKmer length for index hash look up  (default: 32)\n" +
			"  -w, --window\t\tOffset sliding window step size for hash set up (default 32)\n" +
			"  -t, --threads\t\tNumber of CPUs for Goroutine concurrency (default: 4)\n" +
			"  -l, --liftover\tConvert alignment coordinate to linear reference in sam format\n" +
			"  -m, --matrix\t\tScores used to align matches and mismatches (default: humanChimp)\n\n")
			//"Matrix options:\n\n" +
			//"humanChimp:\n" + printMatrix(align.HumanChimpTwoScoreMatrix) + "\n" +
			//"hoxD55:\n" + printMatrix(align.HoxD55ScoreMatrix) + "\n" +
			//"mouseRat:\n" + printMatrix(align.MouseRatScoreMatrix) + "\n")
}

func ggtoolsExtend() {
	fmt.Print(
		"Usage:\n" +
			"  gsw ggtools [options] ref\n\n" +
			"Options:\n" +
			"  -v, --vcf\t\tProvide a VCF to create a graph reference (.gg) used in gsw align\n" +
			"  -a, --axt\t\tUse axt generated from UCSC kentUtils to create a VCF\n\n")
}

func viewExtend() {
	fmt.Print(
		"Usage:\n" +
			"  gsw view [options] ref\n\n" +
			"Options:\n" +
			"  -g, --giraf\t\tProvide a giraf alignment with its reference to visualize the alignment\n" +
			"  -s, --sam\t\tProvide a sam alignment with its reference to visualize the alignment\n\n")
}

type ViewExe struct {
	Cmd *flag.FlagSet
	GirafFile string
	SamFile string
}

func ViewArgs() *ViewExe {
	view := &ViewExe{Cmd: flag.NewFlagSet("view", flag.ExitOnError)}
	view.Cmd.StringVar(&view.GirafFile, "giraf", "", "Visualize sequences aligned to graphs in giraf format\n")
	view.Cmd.StringVar(&view.SamFile, "sam", "", "Visualize sequences aligned to graphs in sam format\n\n")
	return view
}

func Init(fs *flag.FlagSet, args []string) {
	fs.Parse(args)
}

func RunViewExe() error {
	view := ViewArgs()
	Init(view.Cmd, os.Args[2:])

	tail := view.Cmd.Args()
	if len(tail) > 1 || !strings.HasSuffix(tail[0], ".gg") {
		flag.PrintDefaults()
		return fmt.Errorf("Error: Apologies, your command prompt was not recognized...\n\n-xoxo GG\n")
	}
	viewAlignmentStdOut(tail[0], view.GirafFile, view.SamFile)
	return nil
}

func viewAlignmentStdOut(ref string, samfile string, gF string) {
	if !strings.HasSuffix(ref, ".gg") {
		log.Fatalf("Error: Apologies, your command prompt was not recognized...\n\n-xoxo GG\n")
	}
	gg := simpleGraph.Read(ref)
	if strings.Compare(samfile, "sam") == 0 {
		samfile, _ := sam.Read(samfile)
		for _, samline := range samfile.Aln {
			log.Printf("%s\n", simpleGraph.ViewGraphAlignment(samline, gg))
		}
	}
}






