// Command Group: "Sequence Evolution & Reconstruction"

// Reconstruct ancient sequences using extant genomes and a newick tree with branch lengths
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

// IlsReconstructSeqSettings defines the usage settings for the ilsReconstructSeq command
type IlsReconstructSeqSettings struct {
	PostProbsFiles string
	ReconFiles     string
	ChromSizesFile string
	OutDir         string
	Precision      float32
}

// ilsReconstructSeqUsage defines the usage statement for the ilsReconstructSeq command
// TODO FINISH THIS DESCRIPTION
func IlsReconstructSeqUsage(IlsReconstructSeqFlags *flag.FlagSet) {
	fmt.Print(
		"ilsReconstructSeq accounts for incomplete lineage sorting in ancestral sequence reconstruction, by weighting the input sequence reconstructions by the posterior probability at that position of the associated topology (out of four). Memory scales by pDNA size, so if working with large genomes, make sure to split by chromosome. " +
			"This program takes as input a txt file containing the addresses of n wig files, and a txt file containing the addresses of n corresponding pFastas. " +
			"This program returns a pfasta file containing a sequence representing the weighted average of the input pFastas." +
			"Usage:\n" +
			"ilsReconstructSeq <posterior-probabilities>.txt <recons>.txt outDir\n" +
			"options:\n")
	IlsReconstructSeqFlags.PrintDefaults()
}

// parseIlsReconstructSeqArgs is the main function of the IlsReconstructSeq command. It parses options and runs the ilsReconstructSeq function.
func parseIlsReconstructSeqArgs() {
	expectedNumArgs := 4
	var err error
	IlsReconstructSeqFlags := flag.NewFlagSet("IlsReconstructSeq", flag.ExitOnError)
	var precision *float64 = IlsReconstructSeqFlags.Float64("precision", 0.001, "Specify the precision to use when checking for valid pDNA base (sums to 1).")

	err = IlsReconstructSeqFlags.Parse(os.Args[2:])
	exception.PanicOnErr(err)
	IlsReconstructSeqFlags.Usage = func() { IlsReconstructSeqUsage(IlsReconstructSeqFlags) }

	if len(IlsReconstructSeqFlags.Args()) != expectedNumArgs {
		IlsReconstructSeqFlags.Usage()
		log.Fatalf("Error: expecting %d arguments, but got %d\n",
			expectedNumArgs, len(IlsReconstructSeqFlags.Args()))
	}

	postProbsFiles := IlsReconstructSeqFlags.Arg(0)
	reconFiles := IlsReconstructSeqFlags.Arg(1)
	chromSizesFile := IlsReconstructSeqFlags.Arg(1)
	outDir := IlsReconstructSeqFlags.Arg(3)

	s := IlsReconstructSeqSettings{
		PostProbsFiles: postProbsFiles,
		ReconFiles:     reconFiles,
		ChromSizesFile: chromSizesFile,
		OutDir:         outDir,
		Precision:      float32(*precision),
	}

	IlsReconstructSeq(s)
}

func IlsReconstructSeq(s IlsReconstructSeqSettings) {
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

	readRecons.Close()
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

	readPostProbs.Close()
	postProbs := make([]map[string]wig.Wig, 0)

	for _, filepath := range postProbsFileLines {
		postProbs = append(postProbs, wig.Read(filepath, s.ChromSizesFile, 0))
	}

	out := []pFasta.PFasta{reconstruct.IlsReconstructSeq(postProbs, recons, s.Precision)}

	pFasta.Write(s.OutDir, out)
}

func main() {
	parseIlsReconstructSeqArgs()
}
