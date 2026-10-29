// Command Group: "Sequence Evolution & Reconstruction"

// Takes the weighted sum of FASTA/pFASTA sequences given a wig file for all positions with probability
// of that position being one of the sequencess
package main

import (
	"bufio"
	"flag"
	"fmt"
	"log"
	"os"

	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fasta/pFasta"
	"github.com/vertgenlab/gonomics/reconstruct"
	"github.com/vertgenlab/gonomics/wig"
)

// SumIlsSeqsSettings defines the usage settings for the iSumIlsSeqs command
type SumIlsSeqsSettings struct {
	PostProbsFiles string
	ReconFiles     string
	ChromSizesFile string
	OutDir         string
	Precision      float32
}

// SumIlsSeqsUsage defines the usage statement for the iSumIlsSeqs command
func SumIlsSeqsUsage(SumIlsSeqsFlags *flag.FlagSet) {
	fmt.Print(
		"iSumIlsSeqs performs ancestral sequence reconstruction while accounting for phylogenetic topological uncertainty resulting from incomplete lineage sorting.\n" +
			"Initial ancestral sequence reconstructions are generated for each possible topology. A final, adjusted ancestral sequence estimate is generated as the\n" +
			"average of ancestral estimates across topologies, weighted by the posterior probability of a topology at each position.\n" +
			"Memory scales by pDNA size, so if working with large genomes, make sure to split by chromosome.\n" +
			"This program takes as input a txt file containing the addresses of n wig files, and a txt file containing the addresses of n corresponding pFastas.\n" +
			"This program returns a pFasta file containing a sequence representing the weighted average of the input pFastas.\n" +
			"Usage:\n" +
			"iSumIlsSeqs <posterior-probabilities>.txt <recons>.txt outDir\n" +
			"options:\n")
	SumIlsSeqsFlags.PrintDefaults()
}

// parseSumIlsSeqsArgs is the main function of the SumIlsSeqs command. It parses options and runs the iSumIlsSeqs function.
func parseSumIlsSeqsArgs() {
	expectedNumArgs := 4
	var err error
	SumIlsSeqsFlags := flag.NewFlagSet("SumIlsSeqs", flag.ExitOnError)
	var precision *float64 = SumIlsSeqsFlags.Float64("precision", 0.001, "Specify the precision to use when checking for valid pDNA base (sums to 1).")
	SumIlsSeqsFlags.Usage = func() { SumIlsSeqsUsage(SumIlsSeqsFlags) }
	err = SumIlsSeqsFlags.Parse(os.Args[2:])

	exception.PanicOnErr(err)

	if len(SumIlsSeqsFlags.Args()) != expectedNumArgs {
		SumIlsSeqsFlags.Usage()
		log.Fatalf("Error: expecting %d arguments, but got %d\n",
			expectedNumArgs, len(SumIlsSeqsFlags.Args()))
	}

	postProbsFiles := SumIlsSeqsFlags.Arg(0)
	reconFiles := SumIlsSeqsFlags.Arg(1)
	chromSizesFile := SumIlsSeqsFlags.Arg(2)
	outDir := SumIlsSeqsFlags.Arg(3)

	s := SumIlsSeqsSettings{
		PostProbsFiles: postProbsFiles,
		ReconFiles:     reconFiles,
		ChromSizesFile: chromSizesFile,
		OutDir:         outDir,
		Precision:      float32(*precision),
	}

	SumIlsSeqs(s)
}

func SumIlsSeqs(s SumIlsSeqsSettings) {
	var fileScanner *bufio.Scanner
	var err error
	var readPostProbs *os.File
	var readRecons *os.File

	readRecons, err = os.Open(s.ReconFiles)
	if err != nil {
		log.Fatalf("%s does not exist.", s.ReconFiles)
	}

	fileScanner = bufio.NewScanner(readRecons)
	fileScanner.Split(bufio.ScanLines)
	reconFileLines := make([]string, 0)

	for fileScanner.Scan() {
		reconFileLines = append(reconFileLines, fileScanner.Text())
	}

	err = readRecons.Close()
	exception.PanicOnErr(err)
	recons := make([]pFasta.PFasta, 0)

	for _, filepath := range reconFileLines {
		recons = append(recons, pFasta.Read(filepath)[0])
	}

	readPostProbs, err = os.Open(s.PostProbsFiles)
	if err != nil {
		log.Fatalf("%s does not exist.", s.PostProbsFiles)
	}

	fileScanner = bufio.NewScanner(readPostProbs)
	fileScanner.Split(bufio.ScanLines)
	postProbsFileLines := make([]string, 0)

	for fileScanner.Scan() {
		postProbsFileLines = append(postProbsFileLines, fileScanner.Text())
	}

	err = readPostProbs.Close()
	exception.PanicOnErr(err)
	postProbs := make([]map[string]wig.Wig, 0)

	for _, filepath := range postProbsFileLines {
		postProbs = append(postProbs, wig.Read(filepath, s.ChromSizesFile, 0))
	}

	out := []pFasta.PFasta{reconstruct.SumIlsSeqs(postProbs, recons, s.Precision)}

	pFasta.Write(s.OutDir, out)
}

func main() {
	parseSumIlsSeqsArgs()
}
